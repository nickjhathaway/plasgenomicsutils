"""Tests for the Fws within-host diversity core and its two front-ends.

The regression estimator is validated against moimix::getFws out-of-band (see
docs/fws_moimix_reconciliation.md — max |diff| ~5e-5, i.e. output rounding only);
these tests lock in the properties that keep it faithful.
"""

import numpy as np
import pytest

from plasgenomicsutils.lib.fws import AlleleDepths, compute_fws, read_ad_table


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
    d = read_ad_table(str(ad), ["s1", "s2", "s3", "s4"], snps_only=True)
    assert d.ref.shape == (2, 4)                    # snps_only drops the indel row
    assert d.ref[0].tolist() == [8, 5, 0, 0]
    assert d.nonref[0].tolist() == [2, 5, 0, 0]
    assert d.n_multiallelic == 0
    # without snps_only every row with AD is used, so the indel row is included
    d_all = read_ad_table(str(ad), ["s1", "s2", "s3", "s4"])
    assert d_all.n_sites == 3


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
    samples, d_v = read_ad_vcf(str(vcf))

    ad = tmp_path / "ad.tsv"
    ad.write_text("chr1\t10\tA\tT\t8,2\t5,5\t0,9\nchr1\t30\tG\tC\t7,3\t6,4\t7,3\n")
    d_t = read_ad_table(str(ad), samples)

    assert samples == ["s1", "s2", "s3"]
    for field in ("ref", "depth", "sumsq", "pop"):
        assert np.array_equal(getattr(d_v, field), getattr(d_t, field))
    fv, _ = compute_fws(d_v)
    ft, _ = compute_fws(d_t)
    assert np.allclose(fv, ft, equal_nan=True)


# --------------------------------------------------------------------------- #
#  Multiallelic sites                                                         #
# --------------------------------------------------------------------------- #

# One triallelic SNP. Sample AD is (ref, alt1, alt2): a 50:50 mix of the two ALTs, a
# three-way mix, a sample homozygous for alt2, and one homozygous for ref.
TRI_SAMPLES = ["mix_alt1_alt2", "mix_three", "pure_alt2", "pure_ref"]
TRI_AD = [[0, 10, 10], [10, 5, 5], [0, 0, 20], [20, 0, 0]]
TRI_HW = [0.5, 0.625, 0.0, 0.0]            # 1 - Σ q² per sample


def _tri_vcf(tmp_path):
    cells = "\t".join(f"./.:{','.join(map(str, ad))}" for ad in TRI_AD)
    vcf = tmp_path / "tri.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(TRI_SAMPLES) + "\n"
        f"chr1\t10\t.\tA\tT,G\t.\t.\t.\tGT:AD\t{cells}\n")
    return vcf


def _tri_table(tmp_path):
    ad = tmp_path / "tri.tsv"
    ad.write_text("chr1\t10\tA\tT,G\t" + "\t".join(",".join(map(str, a)) for a in TRI_AD) + "\n")
    return ad


def _hw(d):
    return 1 - d.sumsq / (d.depth * d.depth)


def test_biallelic_depths_reduce_to_the_ref_alt_form_exactly():
    # AlleleDepths built from full AD records and from (ref, alt) matrices must agree to
    # the bit at biallelic sites, so nothing about the classic path changed.
    ref, alt = _cohort()
    a = AlleleDepths.from_ref_alt(ref, alt)
    b = AlleleDepths.from_records([np.stack([ref[i], alt[i]], axis=1)
                                   for i in range(ref.shape[0])], ref.shape[1])
    for field in ("ref", "depth", "sumsq", "pop"):
        assert np.array_equal(getattr(a, field), getattr(b, field))
    # and the generalised heterozygosity is the familiar 2q(1-q)
    q = alt / (ref + alt)
    assert np.allclose(_hw(a), 2 * q * (1 - q))
    fa, na = compute_fws(a)
    fb, nb = compute_fws(ref, alt)
    assert np.array_equal(fa, fb) and np.array_equal(na, nb)


def test_ad_table_collapses_a_multiallelic_row_instead_of_truncating_it(tmp_path):
    # the old reader took the first two AD fields of a 'T,G' row, so (0,10,10) became
    # (0,10) -- a perfect 50:50 mixture scored as homozygous. Every allele counts now.
    d = read_ad_table(str(_tri_table(tmp_path)), TRI_SAMPLES)
    assert d.n_sites == 1 and d.n_multiallelic == 1 and d.multiallelic == "collapse"
    assert d.depth[0].tolist() == [20, 20, 20, 20]
    assert d.pop[0].tolist() == [30, 15, 35]
    assert np.allclose(_hw(d)[0], TRI_HW)
    assert "1 multiallelic record(s) collapsed" in d.multiallelic_note()


def test_ad_table_skip_mode_drops_and_counts_a_multiallelic_row(tmp_path):
    d = read_ad_table(str(_tri_table(tmp_path)), TRI_SAMPLES, multiallelic="skip")
    assert d.n_sites == 0 and d.n_multiallelic == 1
    assert "skipped" in d.multiallelic_note()
    with pytest.raises(ValueError, match="multiallelic"):
        read_ad_table(str(_tri_table(tmp_path)), TRI_SAMPLES, multiallelic="first")


def test_snps_only_needs_every_alt_to_be_a_single_base(tmp_path):
    ad = tmp_path / "mixed.tsv"
    ad.write_text("chr1\t10\tA\tT,G\t0,10,10\t10,5,5\n"       # triallelic SNP: kept
                  "chr1\t20\tA\tT,AT\t0,10,10\t10,5,5\n")     # SNP + indel: not a SNP site
    assert read_ad_table(str(ad), ["a", "b"], snps_only=True).n_sites == 1
    assert read_ad_table(str(ad), ["a", "b"], snps_only=False).n_sites == 2


def test_vcf_reader_collapses_multiallelic_records(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.fws import read_ad_vcf

    samples, d = read_ad_vcf(str(_tri_vcf(tmp_path)))
    assert samples == TRI_SAMPLES
    assert d.n_sites == 1 and d.n_multiallelic == 1
    assert np.allclose(_hw(d)[0], TRI_HW)
    # the table front-end reads the same site identically
    d_t = read_ad_table(str(_tri_table(tmp_path)), TRI_SAMPLES)
    for field in ("ref", "depth", "sumsq", "pop"):
        assert np.array_equal(getattr(d, field), getattr(d_t, field))
    _, d_skip = read_ad_vcf(str(_tri_vcf(tmp_path)), multiallelic="skip")
    assert d_skip.n_sites == 0 and d_skip.n_multiallelic == 1


def test_a_mixture_of_two_alt_alleles_lowers_fws_only_when_collapsed():
    # A cohort where sample 0 is monoclonal at every biallelic site but a 50:50 mix of the
    # two ALTs at every triallelic site. Collapsing sees the mixture; skipping the
    # triallelic sites (or splitting them, which is what `bcftools norm -m-` would do)
    # cannot, and calls the sample monoclonal.
    rng = np.random.default_rng(3)
    n_samples, n_bi, n_tri, dp = 8, 300, 100, 40
    records, split_records = [], []
    for i in range(n_bi):
        alt_carriers = rng.integers(1, n_samples)          # varies the population MAF
        ad = np.zeros((n_samples, 2))
        ad[:, 0] = dp
        ad[rng.choice(n_samples, alt_carriers, replace=False)] = [0, dp]
        ad[0] = [dp, 0]                                    # sample 0 always homozygous ref
        records.append(ad)
        split_records.append(ad)
    for i in range(n_tri):
        ad = np.zeros((n_samples, 3))
        which = rng.integers(0, 3, size=n_samples)         # each sample homozygous for one allele
        ad[np.arange(n_samples), which] = dp
        ad[0] = [0, dp / 2, dp / 2]                        # sample 0: mixed between the two ALTs
        records.append(ad)
        split_records.append(ad[:, [0, 1]])                # what `norm -m-` leaves behind
        split_records.append(ad[:, [0, 2]])
    collapsed = AlleleDepths.from_records(records, n_samples)
    skipped = AlleleDepths.from_records(records[:n_bi], n_samples)
    split = AlleleDepths.from_records(split_records, n_samples)

    f_collapsed, _ = compute_fws(collapsed)
    f_skipped, _ = compute_fws(skipped)
    f_split, _ = compute_fws(split)
    assert f_skipped[0] == pytest.approx(1.0, abs=1e-9)
    assert f_split[0] == pytest.approx(1.0, abs=1e-9)
    assert f_collapsed[0] < 0.9
    # the other samples are monoclonal everywhere and stay so under every treatment
    assert np.allclose(f_collapsed[1:], 1.0, atol=1e-9)


# --------------------------------------------------------------------------- #
#  Depth-based trimming of cohort-wide alleles                                #
# --------------------------------------------------------------------------- #

# A callset joint-called across a larger cohort carries ALTs no sample here supports.
# Row 1: a SNP whose second ALT (an insertion) has no reads -> a SNP site once trimmed.
# Row 2: a triallelic SNP whose third allele has no reads -> biallelic once trimmed.
# Row 3: no ALT reads at all -> dropped as monomorphic.
COHORT_TSV = ("chr1\t10\tA\tT,ATT\t8,2,0\t5,5,0\t0,9,0\n"
              "chr1\t20\tG\tC,A\t7,3,0\t6,4,0\t7,3,0\n"
              "chr1\t30\tC\tG\t9,0\t9,0\t.\n")
COHORT_SAMPLES = ["s1", "s2", "s3"]


def test_trimming_drops_alts_no_sample_has_reads_for_before_classifying(tmp_path):
    ad = tmp_path / "cohort.tsv"
    ad.write_text(COHORT_TSV)
    d = read_ad_table(str(ad), COHORT_SAMPLES, snps_only=True)
    assert d.n_sites == 2                       # rows 1 and 2 kept, row 3 monomorphic
    assert d.n_multiallelic == 0                # row 2 is biallelic in this cohort
    assert d.n_alt_trimmed == 2 and d.n_monomorphic == 1
    assert d.pop.shape[1] == 2                  # nothing wider than biallelic survived
    assert d.pop.tolist() == [[13, 16], [20, 10]]
    note = d.multiallelic_note()
    assert "2 record(s) lost ALT allele(s)" in note and "1 record(s) dropped" in note
    # untrimmed, the same rows read as an indel site (dropped by snps_only) and a
    # triallelic one, and the monomorphic row stays in
    u = read_ad_table(str(ad), COHORT_SAMPLES, snps_only=True, trim=False)
    assert u.n_sites == 2 and u.n_multiallelic == 1 and u.n_alt_trimmed == 0
    assert u.pop.shape[1] == 3


def test_a_read_less_allele_changes_nothing_in_the_estimate(tmp_path):
    # trimming is bookkeeping only: a zero-depth allele has zero frequency everywhere,
    # so Fws is identical whether it is trimmed or carried along.
    ad = tmp_path / "cohort.tsv"
    ad.write_text(COHORT_TSV)
    t = read_ad_table(str(ad), COHORT_SAMPLES, snps_only=False, multiallelic="collapse")
    u = read_ad_table(str(ad), COHORT_SAMPLES, snps_only=False, multiallelic="collapse",
                      trim=False)
    # same polymorphic sites; the untrimmed one also carries the monomorphic row
    assert t.n_sites == 2 and u.n_sites == 3
    ft, _ = compute_fws(t, min_alt_samples=1)
    fu, _ = compute_fws(u, min_alt_samples=1)   # min_alt_samples removes row 3 there too
    assert np.allclose(ft, fu, equal_nan=True)


def test_vcf_reader_trims_the_same_way(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.fws import read_ad_vcf

    vcf = tmp_path / "cohort.vcf"
    rows = []
    for line in COHORT_TSV.splitlines():
        chrom, pos, ref, alt, *ads = line.split("\t")
        cells = "\t".join("./.:." if a == "." else f"./.:{a}" for a in ads)
        rows.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:AD\t{cells}")
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(COHORT_SAMPLES)
        + "\n" + "\n".join(rows) + "\n")
    samples, d_v = read_ad_vcf(str(vcf), snps_only=True)
    ad = tmp_path / "cohort.tsv"
    ad.write_text(COHORT_TSV)
    d_t = read_ad_table(str(ad), COHORT_SAMPLES, snps_only=True)
    assert samples == COHORT_SAMPLES
    for field in ("ref", "depth", "sumsq", "pop"):
        assert np.array_equal(getattr(d_v, field), getattr(d_t, field))
    assert (d_v.n_multiallelic, d_v.n_alt_trimmed, d_v.n_monomorphic) == (0, 2, 1)
