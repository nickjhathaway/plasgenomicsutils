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


# --------------------------------------------------------------------------- #
#  Microhaplotypes: long-format allele tables and the unbinned estimator      #
# --------------------------------------------------------------------------- #

# Three loci, three samples. `mono` carries one haplotype everywhere; `mix` is a 50:50 mix
# of two haplotypes at every locus; `absent` has no row at locus L3 (zero depth there).
ALLELE_TSV = ("library_sample_name\ttarget_name\tseq\treads\n"
              "mono\tL1\tAAA\t100\nmix\tL1\tAAA\t50\nmix\tL1\tAAC\t50\nabsent\tL1\tAAC\t80\n"
              "mono\tL2\tGGG\t120\nmix\tL2\tGGG\t60\nmix\tL2\tGGT\t60\nabsent\tL2\tGGT\t90\n"
              "mono\tL3\tTTT\t100\nmix\tL3\tTTT\t50\nmix\tL3\tTTA\t50\n"
              "mono\tL3\tTTT\t10\n")   # a duplicated row is summed, not a fourth allele


def test_allele_table_reader_builds_one_record_per_locus(tmp_path):
    from plasgenomicsutils.lib.fws import read_allele_table

    t = tmp_path / "alleles.tsv"
    t.write_text(ALLELE_TSV)
    samples, d = read_allele_table(str(t))
    assert samples == ["absent", "mix", "mono"]              # sorted
    assert d.n_sites == 3 and d.n_multiallelic == 3
    assert d.pop.shape == (3, 2)                             # two haplotypes per locus
    assert d.depth[2].tolist() == [0, 100, 110]              # L3: absent=0, mix=100, mono=100+10
    assert d.pop[2].tolist() == [160, 50]                    # TTT 110+50, TTA 50
    hw = 1 - d.sumsq / np.where(d.depth > 0, d.depth, np.nan) ** 2
    assert np.allclose(hw[:, 1], 0.5) and np.allclose(hw[:, 2], 0.0)


def test_allele_table_columns_can_be_renamed(tmp_path):
    from plasgenomicsutils.lib.fws import read_allele_table

    t = tmp_path / "alleles.tsv"
    t.write_text(ALLELE_TSV.replace("library_sample_name", "s").replace("target_name", "loc")
                 .replace("seq", "hap").replace("reads", "n"))
    samples, d = read_allele_table(str(t), sample_col="s", locus_col="loc", allele_col="hap",
                                   reads_col="n")
    assert samples == ["absent", "mix", "mono"] and d.n_sites == 3


def test_unbinned_estimator_is_the_per_site_regression(tmp_path):
    from plasgenomicsutils.lib.fws import read_allele_table

    t = tmp_path / "alleles.tsv"
    t.write_text(ALLELE_TSV)
    samples, d = read_allele_table(str(t))
    fws, n = compute_fws(d, n_bins=0)
    assert fws[samples.index("mono")] == pytest.approx(1.0)
    assert n[samples.index("absent")] == 2                   # no depth at L3
    # by hand: Fws = 1 - Σ Hs·Hw / Σ Hs² over the sample's usable loci
    tot = d.pop.sum(1); Hs = 1 - ((d.pop / tot[:, None]) ** 2).sum(1)
    s = samples.index("mix")
    Hw = 1 - d.sumsq[:, s] / d.depth[:, s] ** 2
    assert fws[s] == pytest.approx(1 - (Hs * Hw).sum() / (Hs * Hs).sum())
    assert fws[s] < 0.6
    # the ratio flavour is 1 - Σ Hw / Σ Hs
    fr, _ = compute_fws(d, estimator="ratio", n_bins=0)
    assert fr[s] == pytest.approx(1 - Hw.sum() / Hs.sum())


def test_unbinned_equals_binned_when_population_het_is_constant():
    # with one Hs everywhere the through-origin slope collapses to mean(Hw)/Hs whether the
    # points are sites or bins, so n_bins=0 and the moimix binning must agree exactly there;
    # elsewhere they are different estimators and are documented as such.
    rng = np.random.default_rng(2)
    n = 40
    alt = np.full((n, 8), 3.0)
    a0 = rng.integers(0, 4, size=n).astype(float)
    alt[:, 0] = a0
    alt[:, 1] = 6.0 - a0
    ref = 10.0 - alt
    assert np.allclose(alt.sum(axis=1) / (ref + alt).sum(axis=1), 0.30)
    binned, _ = compute_fws(ref, alt, n_bins=10)
    unbinned, _ = compute_fws(ref, alt, n_bins=0)
    assert np.allclose(binned, unbinned, atol=1e-12)


# --------------------------------------------------------------------------- #
#  Population frequencies supplied from outside                              #
# --------------------------------------------------------------------------- #


def test_supplying_the_cohorts_own_frequencies_changes_nothing(tmp_path):
    from plasgenomicsutils.lib.fws import read_allele_table, read_pop_freqs, write_pop_freqs

    t = tmp_path / "alleles.tsv"
    t.write_text(ALLELE_TSV)
    samples, d = read_allele_table(str(t))
    own = d.population_freqs()
    assert own["L3"] == {"TTT": pytest.approx(160 / 210), "TTA": pytest.approx(50 / 210)}
    f = tmp_path / "freqs.tsv"
    write_pop_freqs(own, str(f))
    back = read_pop_freqs(str(f))
    for n_bins in (0, 10):
        a, _ = compute_fws(d, n_bins=n_bins)
        b, _ = compute_fws(d, n_bins=n_bins, pop_freqs=back)
        assert np.allclose(a, b, equal_nan=True, atol=1e-9)
    assert compute_fws.last_pop_freq_misses == 0


def test_external_frequencies_drive_hs_and_unlisted_loci_are_dropped(tmp_path):
    from plasgenomicsutils.lib.fws import read_allele_table

    t = tmp_path / "alleles.tsv"
    t.write_text(ALLELE_TSV)
    samples, d = read_allele_table(str(t))
    s = samples.index("mix")
    # a reference population where L1 is nearly fixed, L2 is very diverse (one allele of
    # ours plus two never seen here), and L3 is not listed at all
    ref = {"L1": {"AAA": 0.98, "AAC": 0.02},
           "L2": {"GGG": 0.25, "GGT": 0.25, "GGA": 0.25, "GGC": 0.25},
           }
    f, n = compute_fws(d, n_bins=0, pop_freqs=ref)
    assert compute_fws.last_pop_freq_misses == 1
    assert n[s] == 2                                       # L3 dropped for everyone
    Hs = np.array([1 - (0.98**2 + 0.02**2), 1 - 4 * 0.25**2])
    Hw = np.array([0.5, 0.5])                              # mix is 50:50 at both
    assert f[s] == pytest.approx(1 - (Hs * Hw).sum() / (Hs * Hs).sum())
    # frequencies given as counts are renormalised
    ref_counts = {"L1": {"AAA": 98, "AAC": 2}, "L2": {"GGG": 1, "GGT": 1, "GGA": 1, "GGC": 1}}
    f2, _ = compute_fws(d, n_bins=0, pop_freqs=ref_counts)
    assert f2[s] == pytest.approx(f[s])


def test_vcf_sites_are_keyed_chrom_pos_for_frequencies(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.fws import read_ad_vcf

    samples, d = read_ad_vcf(str(_tri_vcf(tmp_path)))
    assert d.sites == ["chr1:10"] and d.alleles == [["A", "T", "G"]]
    own = d.population_freqs()["chr1:10"]
    assert own == {"A": pytest.approx(30 / 80), "T": pytest.approx(15 / 80), "G": pytest.approx(35 / 80)}
    # a population where the site is biallelic A/T: G carries frequency 0 there, and the
    # G-only sample is scored against that population rather than its own cohort's
    f, n = compute_fws(d, n_bins=0, pop_freqs={"chr1:10": {"A": 0.5, "T": 0.5}})
    assert n.tolist() == [1, 1, 1, 1]
    assert f[samples.index("pure_alt2")] == pytest.approx(1.0)      # homozygous: Hw = 0
    assert f[samples.index("mix_alt1_alt2")] == pytest.approx(1 - 0.5 / 0.5)


# --- k-allele sites and the MAF grid ---------------------------------------------------
#
# moimix's grid is ten bins over [0, 0.5], which is the whole range of a biallelic minor-
# allele fraction. A site with k alleles has `1 - max(p)` up to 1 - 1/k, and used to fall
# off the top: the regression estimator gave every such site one shared overflow bin, the
# ratio estimator clipped them into the 0.45-0.5 bin beside biallelic sites of half their
# heterozygosity. The grid now extends upward in the same width -- and the biallelic part of
# it does not move.


def _with_k_allele_sites(n_multi=40, seed=3):
    rng = np.random.default_rng(seed)
    ref, alt = _cohort()
    n_bi, n_samples = ref.shape
    recs = [np.stack([ref[i], alt[i]], axis=1) for i in range(n_bi)]
    for _ in range(n_multi):
        # four alleles at roughly even depth: 1 - max(p) is about 0.7, well above 0.5
        recs.append(rng.integers(20, 30, size=(n_samples, 4)).astype(float))
    return AlleleDepths.from_records(recs, n_samples), n_bi


def test_a_biallelic_panel_is_binned_exactly_as_before():
    # nothing above 0.5, so no extra bins, so the same edges and the same numbers
    ref, alt = _cohort()
    f, n = compute_fws(ref, alt, estimator="regression")
    f2, n2 = compute_fws(AlleleDepths.from_ref_alt(ref, alt), estimator="regression")
    assert np.array_equal(f, f2) and np.array_equal(n, n2)


def test_k_allele_sites_get_bins_of_their_own_rather_than_the_top_biallelic_one():
    d, n_bi = _with_k_allele_sites()
    for est in ("regression", "ratio"):
        f, n = compute_fws(d, estimator=est, min_alt_samples=1)
        assert np.isfinite(f).all()
        assert (n == d.n_sites).all()          # every site, k-allele ones included, counted
    # the mechanism: with the k-allele sites removed, the biallelic sites alone give the
    # same per-sample regression inputs, i.e. adding sites above 0.5 did not move the
    # biallelic bins. Check by scoring the biallelic subset and comparing it to a run on the
    # full panel from which the k-allele sites contribute only their own (new) bins.
    ref, alt = _cohort()
    f_bi, _ = compute_fws(ref, alt, estimator="regression")
    f_all, _ = compute_fws(d, estimator="regression", min_alt_samples=1)
    # not identical -- the extra bins add points to the regression -- but a mixture that is
    # 40 k-allele sites in 340 cannot swing a monoclonal-vs-polyclonal reading; the
    # correlation across samples has to stay near one
    assert np.corrcoef(f_bi, f_all)[0, 1] > 0.9
