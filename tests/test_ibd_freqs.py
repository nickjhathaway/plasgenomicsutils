"""Allele-frequency counting (global + per-region), including multiallelic + missing."""

import math

import pytest

pytest.importorskip("cyvcf2")
from plasgenomicsutils.lib.ibd_freqs import compute_allele_freqs

_VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n"
    "chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t1/1\n"      # alt counts 0,1,2
    "chr1\t20\t.\tA\tT,C\t.\t.\t.\tGT\t0/1\t1/2\t./.\n"     # multiallelic; s3 missing
)


def _af(df, snp):
    return df.set_index("snp_id")["af"].to_dict()[snp]


def test_global_and_region_allele_freqs(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    s2r = {"s1": "A", "s2": "A", "s3": "B"}
    g, r = compute_allele_freqs(str(vcf), s2r)

    # global: site10 = (0+1+2)/6 = 0.5 ; site20 = (1+2)/4 = 0.75 (s3 missing -> excluded)
    assert _af(g, "chr1:10") == pytest.approx(0.5)
    assert _af(g, "chr1:20") == pytest.approx(0.75)

    ga = r[r.region == "A"]
    gb = r[r.region == "B"]
    assert _af(ga, "chr1:10") == pytest.approx(0.25)   # s1,s2 -> (0+1)/4
    assert _af(ga, "chr1:20") == pytest.approx(0.75)   # (1+2)/4
    assert _af(gb, "chr1:10") == pytest.approx(1.0)    # s3 hom-alt
    assert math.isnan(_af(gb, "chr1:20"))              # s3 missing -> AN 0 -> NaN


def test_zero_based_snp_ids(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    g1, _ = compute_allele_freqs(str(vcf))
    g0, _ = compute_allele_freqs(str(vcf), zero_based=True)
    assert list(g1.snp_id) == ["chr1:10", "chr1:20"]
    assert list(g0.snp_id) == ["chr1:9", "chr1:19"]
