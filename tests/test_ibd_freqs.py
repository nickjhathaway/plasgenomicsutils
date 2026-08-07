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
    assert _af(g, "chr1:9") == pytest.approx(0.5)
    assert _af(g, "chr1:19") == pytest.approx(0.75)

    ga = r[r.group == "A"]
    gb = r[r.group == "B"]
    assert _af(ga, "chr1:9") == pytest.approx(0.25)   # s1,s2 -> (0+1)/4
    assert _af(ga, "chr1:19") == pytest.approx(0.75)   # (1+2)/4
    assert _af(gb, "chr1:9") == pytest.approx(1.0)    # s3 hom-alt
    assert math.isnan(_af(gb, "chr1:19"))              # s3 missing -> AN 0 -> NaN


def test_snp_ids_are_always_zero_based(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    g, _ = compute_allele_freqs(str(vcf))
    # VCF POS 10 and 20 -> the canonical 0-based labels; there is no toggle for 1-based
    assert list(g.snp_id) == ["chr1:9", "chr1:19"]
    assert "pos_vcf" not in g.columns              # opt-in, so the table stays small


def test_with_pos_vcf_adds_the_one_based_column(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    g, _ = compute_allele_freqs(str(vcf), with_pos_vcf=True)
    assert list(g.pos_vcf) == [10, 20]             # what the VCF itself says
    assert list(g.snp_id) == ["chr1:9", "chr1:19"]  # the key is unchanged
