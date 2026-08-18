"""Integration test for the cyvcf2-backed filter_ad_regenotype."""

import subprocess

import pytest

pytest.importorskip("cyvcf2")
from plasgenomicsutils.lib.regenotype import filter_ad_regenotype

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    '##FORMAT=<ID=ADS,Number=1,Type=Integer,Description="summed AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
)


def _query(path, fmt):
    out = subprocess.run(["bcftools", "query", "-f", f"%POS[ {fmt}]\n", path],
                         stdout=subprocess.PIPE, text=True).stdout
    return {r.split()[0]: r.split()[1:] for r in out.splitlines()}


@pytest.mark.skipif(subprocess.run(["bash", "-c", "command -v bcftools"],
                                   stdout=subprocess.DEVNULL).returncode != 0,
                    reason="bcftools not on PATH")
def test_filter_ad_regenotype_cleans_and_regenotypes(tmp_path):
    # S1 100: alt has 1 read (< min_reads 2) -> zeroed -> 0/0, ADS 30->29
    # S2 100: balanced 15,15 -> stays het 0/1
    # S3 100: hom-alt 0,40 -> stays 1/1
    inp = tmp_path / "in.vcf"
    inp.write_text(
        HEADER +
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD:ADS\t0/1:29,1:30\t0/1:15,15:30\t1/1:0,40:40\n")
    out = str(tmp_path / "out.vcf")
    filter_ad_regenotype(str(inp), out, min_reads=2, min_freq=0.01, het_min_af=0.2)

    gt = _query(out, "%GT")
    ad = _query(out, "%AD")
    ads = _query(out, "%ADS")
    assert gt["100"] == ["0/0", "0/1", "1/1"]        # S1 artifact allele dropped -> hom-ref
    assert ad["100"] == ["29,0", "15,15", "0,40"]    # S1 alt zeroed
    assert ads["100"] == ["29", "30", "40"]          # S1 ADS recomputed from cleaned AD


# ---- ADS is derived when the callset does not carry it ----------------------------
# The frequency denominator is FORMAT/ADS. Without it the released version raised
# KeyError on the first record, so a callset straight from `bcftools call` could not be
# re-genotyped without running singleton_filter_add_ads first.

def _no_ads_vcf(tmp_path):
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4"]
    # s2 carries a 15% minor clone the caller called 0/1; s4 is a genuine even mix
    row = ("chr1\t1000\t.\tA\tT\t100\t.\t.\tGT:AD\t"
           "0/0:100,0\t0/1:85,15\t1/1:0,100\t0/1:50,50")
    p = tmp_path / "no_ads.vcf"
    p.write_text("\n".join(hdr + [row]) + "\n")
    return str(p)


def test_ads_is_derived_from_ad_when_absent(tmp_path):
    inp = _no_ads_vcf(tmp_path)
    out = str(tmp_path / "out.vcf")
    filter_ad_regenotype(inp, out, min_freq=0.20)

    from cyvcf2 import VCF
    v = next(iter(VCF(out)))
    assert "ADS" in VCF(out).raw_header
    # ADS is the sum of the cleaned AD, per sample
    assert [int(x) for x in v.format("ADS")[:, 0]] == [100, 85, 100, 100]
    # the 15% clone is gone and the genotype follows it down to hom-ref
    assert [list(a) for a in v.format("AD")] == [[100, 0], [85, 0], [0, 100], [50, 50]]
    assert [g[:2] for g in v.genotypes] == [[0, 0], [0, 0], [1, 1], [0, 1]]


def test_deriving_ads_matches_having_computed_it_upfront(tmp_path):
    """The derived denominator has to be the one singleton_filter_add_ads would write."""
    from cyvcf2 import VCF

    inp = _no_ads_vcf(tmp_path)
    derived = str(tmp_path / "derived.vcf")
    filter_ad_regenotype(inp, derived, min_freq=0.20)

    # the same input with ADS already present
    text = open(inp).read().replace(
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        '##FORMAT=<ID=ADS,Number=1,Type=Integer,Description="ADS">')
    text = text.replace("GT:AD\t0/0:100,0\t0/1:85,15\t1/1:0,100\t0/1:50,50",
                        "GT:AD:ADS\t0/0:100,0:100\t0/1:85,15:100\t1/1:0,100:100\t0/1:50,50:100")
    pre = tmp_path / "pre.vcf"
    pre.write_text(text)
    upfront = str(tmp_path / "upfront.vcf")
    filter_ad_regenotype(str(pre), upfront, min_freq=0.20)

    a, b = next(iter(VCF(derived))), next(iter(VCF(upfront)))
    assert [list(x) for x in a.format("AD")] == [list(x) for x in b.format("AD")]
    assert [g[:2] for g in a.genotypes] == [g[:2] for g in b.genotypes]


def test_no_add_ads_passes_records_through_instead_of_raising(tmp_path):
    from cyvcf2 import VCF

    inp = _no_ads_vcf(tmp_path)
    out = str(tmp_path / "out.vcf")
    filter_ad_regenotype(inp, out, min_freq=0.20, add_ads=False)   # must not raise

    v = next(iter(VCF(out)))
    assert [list(a) for a in v.format("AD")] == [[100, 0], [85, 15], [0, 100], [50, 50]]
    assert [g[:2] for g in v.genotypes] == [[0, 0], [0, 1], [1, 1], [0, 1]]
