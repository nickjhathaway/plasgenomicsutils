"""vcf_to_bed: VCF POS is 1-based, BED is not, and that is the entire job.

The test that matters is the round trip -- feed the BED back to `bcftools -R` and it must
select exactly the records it was made from. An off-by-one passes every column-shape check
and fails that one.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import vcf_filters as F

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"),
                                reason="bcftools not on PATH")


def _vcf(tmp_path):
    """A SNP, a 4-base deletion, an MNP and another SNP."""
    rows = [(1000, "A", "T"), (2000, "ATTT", "A"), (3000, "AC", "GT"), (4000, "G", "C")]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"]
    for pos, ref, alt in rows:
        hdr.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t222\t.\t.\tGT\t0/1")
    p = tmp_path / "in.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _rows(path):
    return [ln.split("\t") for ln in open(path).read().splitlines() if ln]


def test_positions_are_zero_based_and_span_the_ref_allele(tmp_path):
    out = str(tmp_path / "o.bed")
    F.vcf_to_bed(_vcf(tmp_path), out)
    r = _rows(out)
    assert [x[1] for x in r] == ["999", "1999", "2999", "3999"]     # POS - 1
    # a SNP is one base; the deletion spans its 4-base REF; the MNP its 2
    assert [int(x[2]) - int(x[1]) for x in r] == [1, 4, 2, 1]


def test_the_bed_selects_exactly_the_records_it_came_from(tmp_path):
    """The round trip. A BED that is one base out still looks perfectly well formed."""
    plain = _vcf(tmp_path)
    vcf = str(tmp_path / "in.vcf.gz")          # -R reads through an index, so make one
    subprocess.run(["bcftools", "view", plain, "-Oz", "-o", vcf], check=True,
                   stderr=subprocess.DEVNULL)
    subprocess.run(["bcftools", "index", "-f", vcf], check=True, stderr=subprocess.DEVNULL)
    bed = str(tmp_path / "o.bed")
    F.vcf_to_bed(vcf, bed, name_column=False)

    def positions(*args):
        out = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\n", *args, vcf],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        return sorted(out.stdout.split("\n"))

    assert positions("-R", bed) == positions()

    # the same check against a deliberately shifted BED, so the assertion above is known to
    # discriminate rather than to pass whatever it is given
    shifted = tmp_path / "shifted.bed"
    shifted.write_text("".join(f"{c}\t{int(s) + 1}\t{int(e) + 1}\n" for c, s, e in _rows(bed)))
    assert positions("-R", str(shifted)) != positions()


def test_snps_only_keeps_substitutions(tmp_path):
    out = str(tmp_path / "o.bed")
    F.vcf_to_bed(_vcf(tmp_path), out, snps_only=True)
    assert [x[1] for x in _rows(out)] == ["999", "3999"]


def test_the_name_column_is_the_canonical_snp_label(tmp_path):
    from plasgenomicsutils.lib.intervals import snp_label

    out = str(tmp_path / "o.bed")
    F.vcf_to_bed(_vcf(tmp_path), out)
    r = _rows(out)
    assert [x[3] for x in r][0] == snp_label("chr1", 999)
    # ...and can be left off for a bare 3-column BED
    bare = str(tmp_path / "b.bed")
    F.vcf_to_bed(_vcf(tmp_path), bare, name_column=False)
    assert all(len(x) == 3 for x in _rows(bare))


def test_no_out_path_writes_to_stdout(tmp_path, capfd):
    F.vcf_to_bed(_vcf(tmp_path))
    out = capfd.readouterr().out
    assert out.startswith("chr1\t999\t1000\t")
    assert len(out.strip().splitlines()) == 4


def test_snp_bed_is_this_restricted_to_snps(tmp_path):
    """The IBD panel writer delegates here, so the two cannot drift apart."""
    a, b = str(tmp_path / "a.bed"), str(tmp_path / "b.bed")
    F.snp_bed(_vcf(tmp_path), a)
    F.vcf_to_bed(_vcf(tmp_path), b, snps_only=True)
    assert open(a).read() == open(b).read()
