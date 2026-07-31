"""Tests for the region filters and the biallelic-SNP (trim) step.

These shell out to bcftools/bedtools; they skip cleanly where those aren't installed.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import vcf_filters as F
from plasgenomicsutils.lib import assets

HAVE_BCFTOOLS = shutil.which("bcftools") is not None
HAVE_BEDTOOLS = shutil.which("bedtools") is not None
needs_bcftools = pytest.mark.skipif(not HAVE_BCFTOOLS, reason="bcftools not on PATH")
needs_both = pytest.mark.skipif(not (HAVE_BCFTOOLS and HAVE_BEDTOOLS),
                                reason="bcftools+bedtools not on PATH")

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)


def _vcf(path, body):
    path.write_text(HEADER + body)
    return str(path)


def _count(path):
    out = subprocess.run(["bcftools", "view", "-H", path], stdout=subprocess.PIPE, text=True)
    return len([x for x in out.stdout.splitlines() if x])


def _positions(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\n", path], stdout=subprocess.PIPE, text=True)
    return sorted(int(x) for x in out.stdout.split())


@needs_bcftools
def test_biallelic_snp_filter_trims_artifact_allele(tmp_path):
    # 100: clean biallelic SNP (kept)
    # 200: "multiallelic" A>T,G but G is in NO genotype -> trimmed -> biallelic -> kept
    # 300: genuinely multiallelic (both T and G used) -> stays multiallelic -> dropped
    # 400: indel -> dropped (not a SNP)
    inp = _vcf(tmp_path / "in.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t1/1:0,9\n"
        "chr1\t200\t.\tA\tT,G\t.\t.\t.\tGT:AD\t0/1:5,5,0\t1/1:0,9,0\n"
        "chr1\t300\t.\tA\tT,G\t.\t.\t.\tGT:AD\t0/1:5,5,0\t2/2:0,0,9\n"
        "chr1\t400\t.\tAT\tA\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n")
    out = str(tmp_path / "out.vcf")
    F.biallelic_snp_filter(inp, out, trim=True)
    assert _positions(out) == [100, 200]   # 200 rescued by trim; 300 (multi) & 400 (indel) dropped

    # without trim, 200 stays multiallelic and is dropped
    out2 = str(tmp_path / "out2.vcf")
    F.biallelic_snp_filter(inp, out2, trim=False)
    assert _positions(out2) == [100]


@needs_bcftools
def test_singleton_add_ads_correct_on_multiallelic(tmp_path):
    # ADS must be the sum of ALL AD entries (ref + every ALT), not just ref+alt1.
    inp = _vcf(tmp_path / "in.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:10,6\t1/1:0,9\n"        # biallelic: 16, 9
        "chr1\t200\t.\tA\tT,G\t.\t.\t.\tGT:AD\t1/2:2,7,5\t0/1:4,8,0\n")  # triallelic: 14, 12
    out = str(tmp_path / "out.vcf")
    F.singleton_add_ads(inp, out, min_samples=0)  # keep both sites
    rows = subprocess.run(["bcftools", "query", "-f", "%POS[ %ADS]\n", out],
                          stdout=subprocess.PIPE, text=True).stdout.strip().splitlines()
    ads = {r.split()[0]: r.split()[1:] for r in rows}
    assert ads["100"] == ["16", "9"]
    assert ads["200"] == ["14", "12"]   # ref+alt1+alt2, not the old ref+alt1 (would be 9 and 12)


@needs_both
def test_region_filter_keep_and_exclude(tmp_path):
    inp = _vcf(tmp_path / "in.vcf",
        "chr1\t150\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"
        "chr1\t550\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,9\t0/0:9,0\n"
        "chr1\t950\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n")
    bed = tmp_path / "region.bed"
    bed.write_text("chr1\t100\t200\nchr1\t900\t1000\n")  # covers 150 and 950, not 550

    keep = str(tmp_path / "keep.vcf")
    F.region_filter(inp, keep, bed=str(bed), exclude=False)
    assert _positions(keep) == [150, 950]

    drop = str(tmp_path / "drop.vcf")
    F.region_filter(inp, drop, bed=str(bed), exclude=True)
    assert _positions(drop) == [550]


def test_bundled_assets_exist():
    for name in ("pf3d7_core_regions", "pf3d7_paralog_genes", "pf3d7_tandem_repeats"):
        assert name in assets.available_assets()
        p = assets.asset_path("builtin:" + name)
        with open(p) as fh:
            assert fh.readline()  # non-empty
    assert assets.resolve_bed("/some/plain/path.bed") == "/some/plain/path.bed"
    with pytest.raises(SystemExit):
        assets.asset_path("builtin:does_not_exist")
