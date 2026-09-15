"""split_by_meta: partition a callset by a metadata column, then make each piece's tags right.

The point of the command over a bare ``bcftools +split`` is what happens to each output
after the samples are subset: AC/AN/AF/MAF recomputed on the subset, an optional per-group
MAF floor judged on the group's own frequencies, and an optional ALT trim that is off by
default so the pieces stay merge-able.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib.vcf_filters import split_by_meta

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"),
                                reason="bcftools not on PATH")


def _write_vcf(path, samples, rows):
    """rows: list of (alt, [gt-per-sample]). One contig, GT-only."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    lines = list(hdr)
    for i, (alt, gts) in enumerate(rows, start=1):
        cells = "\t".join(f"{g}" for g in gts)
        lines.append(f"chr1\t{i*100}\t.\tA\t{alt}\t50\tPASS\t.\tGT\t{cells}")
    path.write_text("\n".join(lines) + "\n")


def _meta(path, mapping):
    path.write_text("sample\tcountry\n" + "".join(f"{s}\t{c}\n" for s, c in mapping.items()))


def _samples(path):
    return subprocess.run(["bcftools", "query", "-l", str(path)],
                          capture_output=True, text=True).stdout.split()


def _af(path):
    """{(pos): AF-string} straight from the record, to prove tags were refilled."""
    out = subprocess.run(["bcftools", "query", "-f", "%POS\t%INFO/AF\n", str(path)],
                         capture_output=True, text=True).stdout
    return dict(l.split("\t") for l in out.splitlines())


def test_splits_by_group_and_refills_tags(tmp_path):
    samples = ["g1a", "g1b", "g2a", "g2b"]
    # site 1: variant only in Ghana; site 2: variant only in Kenya; site 3: shared
    rows = [("T", ["0/1", "1/1", "0/0", "0/0"]),
            ("C", ["0/0", "0/0", "0/1", "1/1"]),
            ("G", ["0/1", "0/0", "0/1", "0/0"])]
    vcf = tmp_path / "in.vcf"
    meta = tmp_path / "meta.tsv"
    _write_vcf(vcf, samples, rows)
    _meta(meta, {"g1a": "Ghana", "g1b": "Ghana", "g2a": "Kenya", "g2b": "Kenya"})

    outdir = tmp_path / "out"
    paths = split_by_meta(str(vcf), str(outdir), meta=str(meta), group_col="country")

    assert set(paths) == {"Ghana", "Kenya"}
    assert _samples(paths["Ghana"]) == ["g1a", "g1b"]
    assert _samples(paths["Kenya"]) == ["g2a", "g2b"]
    # Ghana site 1: 3 ALT copies of 4 alleles -> AF 0.75 on the SUBSET, not the parent's 0.375
    assert _af(paths["Ghana"])["100"] == "0.75"


def test_per_group_maf_floor_is_judged_within_the_group(tmp_path):
    samples = ["g1a", "g1b", "g2a", "g2b"]
    # site at pos 100 is polymorphic in Ghana (AF .25) but monomorphic-ref in Kenya
    rows = [("T", ["0/1", "0/0", "0/0", "0/0"])]
    vcf = tmp_path / "in.vcf"
    meta = tmp_path / "meta.tsv"
    _write_vcf(vcf, samples, rows)
    _meta(meta, {"g1a": "Ghana", "g1b": "Ghana", "g2a": "Kenya", "g2b": "Kenya"})

    outdir = tmp_path / "out"
    paths = split_by_meta(str(vcf), str(outdir), meta=str(meta), group_col="country",
                          maf_min=0.1)
    # kept for Ghana (MAF .25 >= .1), dropped for Kenya (monomorphic there)
    n = lambda p: len(subprocess.run(["bcftools", "view", "-H", str(p)],
                                     capture_output=True, text=True).stdout.splitlines())
    assert n(paths["Ghana"]) == 1
    assert n(paths["Kenya"]) == 0


def test_trim_alts_is_optional_and_prunes_group_unused_alleles(tmp_path):
    samples = ["g1a", "g1b", "g2a", "g2b"]
    # multiallelic: ALT allele 2 (C) is used only by Kenya; Ghana uses only allele 1 (T)
    rows = [("T,C", ["0/1", "1/1", "0/2", "2/2"])]
    vcf = tmp_path / "in.vcf"
    meta = tmp_path / "meta.tsv"
    _write_vcf(vcf, samples, rows)
    _meta(meta, {"g1a": "Ghana", "g1b": "Ghana", "g2a": "Kenya", "g2b": "Kenya"})

    def alts(p):
        return subprocess.run(["bcftools", "query", "-f", "%ALT\n", str(p)],
                              capture_output=True, text=True).stdout.strip()

    # default: full ALT set retained (merge-able)
    kept = split_by_meta(str(vcf), str(tmp_path / "keep"), meta=str(meta),
                         group_col="country")
    assert alts(kept["Ghana"]) == "T,C"

    # trim: Ghana keeps only T, Kenya keeps only C
    trimmed = split_by_meta(str(vcf), str(tmp_path / "trim"), meta=str(meta),
                            group_col="country", trim_alts=True)
    assert alts(trimmed["Ghana"]) == "T"
    assert alts(trimmed["Kenya"]) == "C"
