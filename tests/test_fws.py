"""Tests for the Fws within-host diversity core and its two front-ends.

The regression estimator is validated against moimix::getFws out-of-band (see
docs/fws_moimix_reconciliation.md — max |diff| ~5e-5, i.e. output rounding only);
these tests lock in the properties that keep it faithful.
"""

import numpy as np
import pytest

from plasgenomicsutils.lib.fws import compute_fws, read_ad_table


def _cohort():
    # a reproducible cohort whose population MAF spans several bins
    rng = np.random.default_rng(1)
    n = 300
    ref = rng.integers(3, 60, size=(n, 6)).astype(float)
    alt = rng.integers(3, 60, size=(n, 6)).astype(float)
    return ref, alt


def test_monoclonal_sample_scores_one():
    # a pure-homozygous sample (all reads on the reference allele) has zero within-sample
    # heterozygosity everywhere, so the regression slope is 0 and Fws == 1.
    ref, alt = _cohort()
    ref[:, 0] = 30.0
    alt[:, 0] = 0.0
    fws, n_info = compute_fws(ref, alt, estimator="regression")
    assert fws[0] == pytest.approx(1.0, abs=1e-9)
    assert n_info[0] > 0
    fws_r, _ = compute_fws(ref, alt, estimator="ratio", min_alt_samples=2)
    assert fws_r[0] == pytest.approx(1.0, abs=1e-9)


def test_regression_and_ratio_estimators_differ():
    # the two estimators are genuinely different (regression weights bins by pop-het^2),
    # so the same data yields different Fws — a guard against silently swapping them.
    ref, alt = _cohort()
    reg, _ = compute_fws(ref, alt, estimator="regression")
    rat, _ = compute_fws(ref, alt, estimator="ratio")
    assert np.nanmax(np.abs(reg - rat)) > 1e-3


def test_regression_single_bin_reduces_to_ratio_of_means():
    # when every site shares one MAF bin, the through-origin slope Σxy/Σx² collapses to
    # y/x = mean(Hw)/mean(Hs); confirm regression == the hand ratio for that degenerate case.
    rng = np.random.default_rng(2)
    n = 40
    # Hold the population alt total constant across sites (every site pop freq == 0.30, one
    # bin) while sample 0 still varies within-sample: sample 1 compensates sample 0's shift,
    # so Σ alt is fixed. 8 samples x depth 10; sites 2..7 fixed at 3 alt.
    alt = np.full((n, 8), 3.0)
    a0 = rng.integers(0, 4, size=n).astype(float)   # sample 0 alt in 0..3
    alt[:, 0] = a0
    alt[:, 1] = 6.0 - a0                             # so alt0 + alt1 == 6 (constant)
    ref = 10.0 - alt                                 # every sample depth 10
    depth = ref + alt
    p_site = alt.sum(axis=1) / depth.sum(axis=1)     # per-site population alt freq
    assert np.allclose(p_site, 0.30)                 # guard: truly a single MAF bin
    fws, _ = compute_fws(ref, alt, estimator="regression", n_bins=10)
    Hs_site = 2 * p_site * (1 - p_site)              # 0.42 at every site
    q0 = alt[:, 0] / depth[:, 0]
    Hw0 = 2 * q0 * (1 - q0)
    expected = 1 - Hw0.mean() / Hs_site.mean()       # single shared bin -> slope = y/x
    assert fws[0] == pytest.approx(expected, abs=1e-9)


def test_read_ad_table_parses_and_guards(tmp_path):
    ad = tmp_path / "ad.tsv"
    ad.write_text(
        "chr1\t10\tA\tT\t8,2\t5,5\t.\t0,0\n"      # biallelic SNP, one '.' -> (0,0)
        "chr1\t20\tAC\tA\t9,1\t9,1\t9,1\t9,1\n"    # indel REF -> skipped
        "chr1\t30\tG\tC\t7,3\t6,4\t7,3\t6,4\n")
    ref, alt = read_ad_table(str(ad), ["s1", "s2", "s3", "s4"], snps_only=True)
    assert ref.shape == (2, 4)                      # snps_only drops the indel row
    assert ref[0].tolist() == [8, 5, 0, 0]
    assert alt[0].tolist() == [2, 5, 0, 0]
    # default keeps all biallelic rows (moimix parity), so the indel row is included
    ref_all, _ = read_ad_table(str(ad), ["s1", "s2", "s3", "s4"])
    assert ref_all.shape == (3, 4)


def test_read_ad_table_raises_on_total_column_mismatch(tmp_path):
    ad = tmp_path / "bad.tsv"
    ad.write_text("chr1\t10\tA\tT\t8,2\t5,5\n")     # 2 AD columns
    with pytest.raises(ValueError, match="misaligned"):
        read_ad_table(str(ad), ["s1", "s2", "s3"])  # but 3 samples asserted


def test_vcf_and_ad_table_frontends_agree(tmp_path):
    # the two front-ends must yield identical depth matrices from the same data, so a VCF
    # and its bcftools-query AD table give the same Fws.
    cyvcf2 = pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.fws import read_ad_vcf

    vcf = tmp_path / "mini.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n"
        "chr1\t10\t.\tA\tT\t.\t.\t.\tGT:AD\t0/0:8,2\t0/1:5,5\t1/1:0,9\n"
        "chr1\t30\t.\tG\tC\t.\t.\t.\tGT:AD\t0/1:7,3\t0/0:6,4\t0/1:7,3\n")
    samples, ref_v, alt_v = read_ad_vcf(str(vcf))

    ad = tmp_path / "ad.tsv"
    ad.write_text("chr1\t10\tA\tT\t8,2\t5,5\t0,9\nchr1\t30\tG\tC\t7,3\t6,4\t7,3\n")
    ref_t, alt_t = read_ad_table(str(ad), samples)

    assert samples == ["s1", "s2", "s3"]
    assert np.array_equal(ref_v, ref_t)
    assert np.array_equal(alt_v, alt_t)
    fv, _ = compute_fws(ref_v, alt_v)
    ft, _ = compute_fws(ref_t, alt_t)
    assert np.allclose(fv, ft, equal_nan=True)
