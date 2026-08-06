"""Multiple-testing corrections for the IBD selection statistic.

Both a family-wise (Bonferroni) and a false-discovery-rate (Benjamini-Hochberg) cutoff are
reported, along with the genomic inflation factor that says whether either can be believed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import pytest
from scipy.stats import chi2


# --------------------------------------------------------------------------- #
#  Multiple testing                                                            #
# --------------------------------------------------------------------------- #


def test_bh_matches_the_textbook_definition():
    from plasgenomicsutils.lib.ibd_selection import benjamini_hochberg

    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216])
    q = benjamini_hochberg(p)
    m = len(p)
    expected = np.minimum.accumulate((p * m / np.arange(1, m + 1))[::-1])[::-1]
    assert q == pytest.approx(expected)
    assert np.all(np.diff(q) >= -1e-12)          # monotone, by construction
    assert np.all(q >= p)                        # never smaller than the raw p-value


def test_bh_is_less_strict_than_bonferroni_but_ordered_the_same():
    from plasgenomicsutils.lib.ibd_selection import benjamini_hochberg

    rng = np.random.default_rng(0)
    p = np.concatenate([rng.uniform(0, 1, 990), rng.uniform(0, 1e-6, 10)])
    q = benjamini_hochberg(p)
    assert (q < 0.05).sum() >= (p < 0.05 / len(p)).sum()
    # q is monotone in p (it ties in places, so the orderings are not element-identical)
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)


def test_bh_keeps_missing_values_missing():
    from plasgenomicsutils.lib.ibd_selection import benjamini_hochberg

    q = benjamini_hochberg(np.array([0.01, np.nan, 0.5]))
    assert np.isnan(q[1])
    assert np.all(np.isfinite(q[[0, 2]]))
    assert np.all(np.isnan(benjamini_hochberg(np.array([np.nan, np.nan]))))


def test_genomic_inflation_is_one_for_a_true_null():
    from plasgenomicsutils.lib.ibd_selection import genomic_inflation

    rng = np.random.default_rng(1)
    chi2_null = rng.chisquare(1, 200_000)
    assert genomic_inflation(chi2_null) == pytest.approx(1.0, abs=0.02)
    # the real statistic is nothing like this: a tight bulk with heavy tails reads far low
    tight = np.concatenate([rng.chisquare(1, 199_000) * 0.1, rng.chisquare(1, 1000) * 50])
    assert genomic_inflation(tight) < 0.3


def test_assemble_output_reports_both_corrections_and_the_diagnostic():
    from plasgenomicsutils.lib.ibd_selection import assemble_output

    rng = np.random.default_rng(2)
    n = 5000
    z = rng.normal(size=n)
    z[:5] = 9.0                                   # a handful of real signals
    stats = {"maf": rng.uniform(0.05, 0.5, n), "bin_id": np.zeros(n, int),
             "raw_stat": z, "z_score": z, "chi2_stat": z**2,
             "pval": chi2.sf(z**2, 1), "neg_log10_p": -np.log10(chi2.sf(z**2, 1))}
    snp = pd.DataFrame({"snp_id": [f"c:{i}" for i in range(n)],
                        "chr": "c", "pos": np.arange(n)})
    out, info = assemble_output(snp, stats, alpha=0.05)

    assert {"q_value", "significant", "significant_fdr"} <= set(out.columns)
    assert info["n_significant"] == 5             # Bonferroni finds the planted signals
    assert info["n_significant_fdr"] >= info["n_significant"]   # FDR is never stricter
    assert info["lambda_gc"] == pytest.approx(1.0, abs=0.1)     # this null *is* calibrated
    # the BH line is the k*q/m critical value, never stricter than Bonferroni's k = 1 case
    assert info["neg_log10_p_fdr_threshold"] <= info["neg_log10_p_threshold"]
    # and the drawn line reproduces the flag it is drawn beside
    assert ((out["neg_log10_p"] >= info["neg_log10_p_fdr_threshold"]).sum()
            == info["n_significant_fdr"])


def test_fdr_alpha_moves_only_the_fdr_flag():
    from plasgenomicsutils.lib.ibd_selection import assemble_output

    rng = np.random.default_rng(3)
    n = 2000
    z = rng.normal(size=n); z[:20] = 6.0
    stats = {"maf": rng.uniform(0.05, 0.5, n), "bin_id": np.zeros(n, int),
             "raw_stat": z, "z_score": z, "chi2_stat": z**2,
             "pval": chi2.sf(z**2, 1), "neg_log10_p": -np.log10(chi2.sf(z**2, 1))}
    snp = pd.DataFrame({"snp_id": [f"c:{i}" for i in range(n)],
                        "chr": "c", "pos": np.arange(n)})
    strict = assemble_output(snp, stats, alpha=0.05, fdr_alpha=0.01)[1]
    loose = assemble_output(snp, stats, alpha=0.05, fdr_alpha=0.20)[1]
    assert strict["n_significant"] == loose["n_significant"]         # Bonferroni unmoved
    assert strict["n_significant_fdr"] <= loose["n_significant_fdr"]


# --------------------------------------------------------------------------- #
#  Permutation threshold                                                       #
# --------------------------------------------------------------------------- #


def _ibd_matrix(n_pairs=300, n_snps=2000, seg=40, hotspot=None, frac=0.8, seed=0):
    """Pairs with one IBD segment each, scattered at random -- or, with `hotspot`, a
    `frac` share of them all covering the same locus, which is the signal the statistic
    exists to find."""
    from scipy import sparse

    rng = np.random.default_rng(seed)
    rows, cols = [], []
    for i in range(n_pairs):
        start = (hotspot if (hotspot is not None and i < int(frac * n_pairs))
                 else int(rng.integers(0, n_snps - seg)))
        rows.extend([i] * seg)
        cols.extend(range(start, start + seg))
    return sparse.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                             shape=(n_pairs, n_snps)).tocsr()


def _af(n_snps=2000, seed=0):
    """Varied allele frequencies, so the MAF binning is not degenerate."""
    return np.random.default_rng(seed).uniform(0.05, 0.5, n_snps)


def test_the_permutation_keeps_each_pair_s_sharing_and_moves_only_its_position():
    from scipy import sparse

    mat = _ibd_matrix()
    coo = mat.tocoo()
    rng = np.random.default_rng(1)
    shift = rng.integers(0, mat.shape[1], size=mat.shape[0])
    col = (coo.col + shift[coo.row]) % mat.shape[1]
    null = sparse.coo_matrix((coo.data, (coo.row, col)), shape=mat.shape).tocsr()

    # every pair shares exactly as much as before -- relatedness is not what is permuted
    assert np.array_equal(np.asarray(mat.sum(axis=1)).ravel(),
                          np.asarray(null.sum(axis=1)).ravel())
    # ...but the per-SNP totals have moved
    assert not np.array_equal(np.asarray(mat.sum(axis=0)).ravel(),
                              np.asarray(null.sum(axis=0)).ravel())


def _snp_frame(n):
    return pd.DataFrame({"snp_id": [f"c:{i}" for i in range(n)],
                         "chr": "c", "pos": np.arange(n)})


def test_the_threshold_sits_above_a_no_signal_scan_and_below_a_real_one():
    from plasgenomicsutils.lib.ibd_selection import (compute_selection_statistic,
                                                     permutation_null)
    af = _af()
    top = lambda m: np.nanmax(compute_selection_statistic(m, af, n_bins=10)[0]["neg_log10_p"])

    # scattered segments: the observed max is nothing its own null cannot produce, even
    # though plenty of SNPs would clear a Bonferroni line
    flat = _ibd_matrix(hotspot=None, seed=2)
    assert top(flat) < permutation_null(flat, af, n_perm=40, n_bins=10, seed=3)["threshold"]

    # most pairs sharing one locus is the signal, and it clears its own null
    spike = _ibd_matrix(hotspot=900, seed=2)
    assert top(spike) > permutation_null(spike, af, n_perm=40, n_bins=10, seed=3)["threshold"]


def test_the_permutation_threshold_is_reproducible_and_flows_into_the_output():
    from plasgenomicsutils.lib.ibd_selection import (assemble_output,
                                                     compute_selection_statistic,
                                                     permutation_null)
    mat = _ibd_matrix(hotspot=900, seed=4)
    af = _af()
    run = lambda seed: permutation_null(mat, af, n_perm=20, n_bins=10, seed=seed)
    a = run(11)
    assert a["threshold"] == run(11)["threshold"]         # same seed, same threshold
    assert a["threshold"] != run(12)["threshold"]

    stats, _ = compute_selection_statistic(mat, af, n_bins=10)
    out, info = assemble_output(_snp_frame(mat.shape[1]), stats, alpha=0.05, perm=a)
    assert "significant_perm" in out.columns
    assert info["neg_log10_p_perm_threshold"] == pytest.approx(a["threshold"])
    assert info["n_significant_perm"] == int((out["neg_log10_p"] >= a["threshold"]).sum())
    # and it is the strictest of the three on this data
    assert info["n_significant_perm"] <= info["n_significant"] <= info["n_significant_fdr"]


def test_empirical_p_values_are_calibrated_and_ordered_like_the_statistic():
    from plasgenomicsutils.lib.ibd_selection import permutation_null

    mat = _ibd_matrix(hotspot=900, seed=6)
    af = _af()
    n_perm, n_snps = 30, mat.shape[1]
    perm = permutation_null(mat, af, n_perm=n_perm, n_bins=10, seed=7)

    for k in ("p_pointwise", "p_pooled"):
        v = perm[k][np.isfinite(perm[k])]
        assert v.min() > 0 and v.max() <= 1          # Phipson-Smyth: never exactly 0

    # the pooled p-value resolves far finer, which is the whole reason it can feed BH
    assert perm["p_pointwise"][np.isfinite(perm["p_pointwise"])].min() == pytest.approx(
        1 / (1 + n_perm))
    assert perm["p_pooled"][np.isfinite(perm["p_pooled"])].min() < 1 / (1 + n_perm)
    assert perm["n_pool"] == pytest.approx(n_perm * n_snps, rel=0.01)

    # A permutation recalibrates the axis; it does not reorder the SNPs. The empirical
    # p is a step function of the statistic, so equal p for unequal scores is expected --
    # what must never happen is an inversion.
    from plasgenomicsutils.lib.ibd_selection import compute_selection_statistic
    obs = compute_selection_statistic(mat, af, n_bins=10)[0]["neg_log10_p"]
    ok = np.isfinite(obs) & np.isfinite(perm["p_pooled"])
    by_score = np.argsort(obs[ok])
    assert np.all(np.diff(perm["p_pooled"][ok][by_score]) <= 0)   # non-increasing
    # not exactly -1: under the upper tail every deficit shares p ~ 1, so the empirical
    # p-values carry a large tie block that Spearman penalises
    assert spearmanr(obs[ok], perm["p_pooled"][ok]).statistic < -0.99


def test_empirical_fdr_lands_between_the_family_wise_line_and_plain_bh():
    from plasgenomicsutils.lib.ibd_selection import (assemble_output,
                                                     compute_selection_statistic,
                                                     permutation_null)
    mat = _ibd_matrix(hotspot=900, seed=8)
    af = _af()
    perm = permutation_null(mat, af, n_perm=300, n_bins=10, seed=9)
    stats, _ = compute_selection_statistic(mat, af, n_bins=10)
    out, info = assemble_output(_snp_frame(mat.shape[1]), stats, alpha=0.05, perm=perm)

    for c in ("p_pointwise", "p_empirical", "q_empirical", "significant_fdr_perm"):
        assert c in out.columns
    assert info["q_empirical_floor"] < 0.05 / 10          # enough replicates to resolve
    # controlling FDR is looser than controlling the family-wise error, on either null
    assert info["n_significant_perm"] <= info["n_significant_fdr_perm"]
    # ...and the empirical pair is stricter than the chi2(1) pair it replaces
    assert info["n_significant_fdr_perm"] <= info["n_significant_fdr"]
    # the reported line reproduces the flag
    line = info["neg_log10_p_emp_fdr_threshold"]
    assert np.isfinite(line)
    assert int((out["neg_log10_p"] >= line).sum()) == info["n_significant_fdr_perm"]
    # the exchangeability check ran over every bin
    assert np.isfinite(info["perm_bin_tail_min"]) and np.isfinite(info["perm_bin_tail_max"])
    assert info["n_perm"] == 300


def test_the_incremental_counts_match_keeping_every_null_value():
    """`permutation_null` accumulates exceedances with a running searchsorted rather than
    storing the whole (n_perm x n_snps) null, which is the one place an off-by-one would
    hide. Replay the same shifts keeping everything, and compare."""
    from scipy import sparse

    from plasgenomicsutils.lib.ibd_selection import (compute_selection_statistic,
                                                     permutation_null)
    mat = _ibd_matrix(n_pairs=80, n_snps=400, seg=20, hotspot=180, seed=3)
    af = _af(400, seed=3)
    n_perm, seed = 25, 17
    fast = permutation_null(mat, af, n_perm=n_perm, n_bins=8, seed=seed)

    coo = mat.tocoo()
    n_pairs, n_snps = mat.shape
    rng = np.random.default_rng(seed)
    obs = compute_selection_statistic(mat, af, n_bins=8)[0]["neg_log10_p"]
    every = np.empty((n_perm, n_snps))
    for r in range(n_perm):
        shift = rng.integers(0, n_snps, size=n_pairs)
        col = (coo.col + shift[coo.row]) % n_snps
        null = sparse.coo_matrix((coo.data, (coo.row, col)), shape=mat.shape).tocsr()
        every[r] = compute_selection_statistic(null, af, n_bins=8)[0]["neg_log10_p"]

    flat = every[np.isfinite(every)]
    ok = np.isfinite(obs)
    point = np.array([np.sum(every[:, j] >= obs[j]) for j in range(n_snps)])
    pool = np.array([np.sum(flat >= obs[j]) if ok[j] else 0 for j in range(n_snps)])

    assert fast["n_pool"] == flat.size
    assert np.allclose(fast["p_pointwise"],
                       np.where(ok, (1 + point) / (1 + n_perm), np.nan), equal_nan=True)
    assert np.allclose(fast["p_pooled"],
                       np.where(ok, (1 + pool) / (1 + flat.size), np.nan), equal_nan=True)
    assert fast["threshold"] == pytest.approx(
        np.quantile(np.nanmax(every, axis=1), 0.95))

    # the stratified count compares a SNP only against nulls from its own MAF bin
    bins = compute_selection_statistic(mat, af, n_bins=8)[0]["bin_id"]
    strat = np.full(n_snps, np.nan)
    for k in np.unique(bins[bins >= 0]):
        idx = np.where(bins == k)[0]
        pool_k = every[:, idx][np.isfinite(every[:, idx])]
        for j in idx:
            if ok[j]:
                strat[j] = (1 + np.sum(pool_k >= obs[j])) / (1 + pool_k.size)
    assert np.allclose(fast["p_stratified"], strat, equal_nan=True)


def test_stratifying_by_maf_bin_trades_resolution_for_dropping_an_assumption():
    from plasgenomicsutils.lib.ibd_selection import (assemble_output,
                                                     compute_selection_statistic,
                                                     permutation_null)
    mat = _ibd_matrix(hotspot=900, seed=8)
    # A real cohort's allele frequencies are k/n for a smallish n, so MAF is heavily tied
    # and the equal-frequency bin edges collapse: the bins come out very unequal, which is
    # what makes the per-SNP resolution under bin pooling uneven too. A uniform `_af()`
    # would give tidy equal bins and hide that entirely.
    rng = np.random.default_rng(21)
    af = rng.choice([0.02, 0.04, 0.06, 0.12, 0.25, 0.45], size=2000,
                    p=[0.45, 0.2, 0.15, 0.1, 0.06, 0.04])
    n_bins = 10
    # few enough replicates that the smaller bins cannot resolve a q of 0.05 while the
    # larger ones can -- the uneven-floor regime a real cohort lands in
    perm = permutation_null(mat, af, n_perm=60, n_bins=n_bins, seed=9)
    stats, _ = compute_selection_statistic(mat, af, n_bins=n_bins)
    snp = _snp_frame(mat.shape[1])

    # both columns are always written; only the q-value follows `pool`
    g_out, g_info = assemble_output(snp, stats, alpha=0.05, perm=perm, pool="global")
    b_out, b_info = assemble_output(snp, stats, alpha=0.05, perm=perm, pool="bin")
    for out in (g_out, b_out):
        assert {"p_empirical", "p_empirical_binned"} <= set(out.columns)
    assert g_out["p_empirical"].equals(b_out["p_empirical"])
    assert not g_out["q_empirical"].equals(b_out["q_empirical"])
    assert g_info["empirical_pool"] == "global" and b_info["empirical_pool"] == "bin"

    # staying in-bin costs resolution, which is the whole reason global pooling is default
    assert b_info["q_empirical_floor"] > g_info["q_empirical_floor"]
    # ...and unevenly: the floor is per-SNP under bin pooling, so a share of SNPs sit in
    # bins too small to reach the cutoff at all. Global pooling gives every SNP the same
    # floor, so that share is all-or-nothing.
    assert b_info["frac_q_unreachable"] > g_info["frac_q_unreachable"]
    assert g_info["frac_q_unreachable"] in (0.0, 1.0)
    assert 0 < b_info["frac_q_unreachable"] < 1

    # the binning is reported as used, not as requested -- MAF ties collapse the edges
    assert 0 < g_info["n_bins_used"] <= n_bins
    assert 0 < g_info["largest_bin_frac"] <= 1

    assemble_output(snp, stats, alpha=0.05, perm=perm, pool="global")   # no error
    with pytest.raises(ValueError, match="pool must be"):
        assemble_output(snp, stats, alpha=0.05, perm=perm, pool="nope")


def test_too_few_replicates_cannot_resolve_an_fdr_and_the_output_says_so():
    """BH multiplies the rank-1 p-value by m, and the empirical p bottoms out at
    1/(1 + n_perm * n_snps) -- so below roughly 1/fdr_alpha replicates the whole column is
    dead. It must report that rather than quietly returning zero discoveries."""
    from plasgenomicsutils.lib.ibd_selection import (assemble_output,
                                                     compute_selection_statistic,
                                                     permutation_null)
    mat = _ibd_matrix(hotspot=900, seed=8)
    af = _af()
    stats, _ = compute_selection_statistic(mat, af, n_bins=10)
    snp = _snp_frame(mat.shape[1])

    thin = permutation_null(mat, af, n_perm=10, n_bins=10, seed=9)
    _, lean = assemble_output(snp, stats, alpha=0.05, perm=thin)
    # a lone top SNP could not clear q<0.05 at this resolution...
    assert lean["q_empirical_floor"] > 0.05
    # ...though a block of SNPs tied at the limit still can, since BH divides by rank --
    # which is exactly why the floor is reported as a bound and not as a verdict
    k = int(np.ceil(lean["q_empirical_floor"] / 0.05))
    assert lean["n_significant_fdr_perm"] == 0 or lean["n_significant_fdr_perm"] >= k
    # the family-wise threshold works at that many replicates either way
    assert lean["n_significant_perm"] >= 1

    # more replicates buy resolution proportionally
    thick = permutation_null(mat, af, n_perm=100, n_bins=10, seed=9)
    _, rich = assemble_output(snp, stats, alpha=0.05, perm=thick)
    assert rich["q_empirical_floor"] == pytest.approx(lean["q_empirical_floor"] / 10, rel=0.05)


def test_no_permutation_leaves_the_columns_out():
    from plasgenomicsutils.lib.ibd_selection import assemble_output, compute_selection_statistic

    mat = _ibd_matrix(seed=5)
    stats, _ = compute_selection_statistic(mat, _af(), n_bins=10)
    out, info = assemble_output(_snp_frame(mat.shape[1]), stats, alpha=0.05)
    for c in ("significant_perm", "p_pointwise", "p_empirical", "q_empirical",
              "significant_fdr_perm"):
        assert c not in out.columns
    assert np.isnan(info["neg_log10_p_perm_threshold"])
    assert info["n_perm"] == 0


# --------------------------------------------------------------------------- #
#  The statistic itself: variant, tail, and the precision guard                #
# --------------------------------------------------------------------------- #


def test_the_statistic_does_not_collapse_when_the_precision_changes():
    """The guard that catches a cancelling formula.

    `variant="published"` centres each SNP and then sums that same SNP, which is exactly
    zero, so what survives is rounding residue: it shrinks by ~1e9 between float32 and
    float64 and the two give unrelated answers. A statistic worth reporting is insensitive
    to the accumulator's width, so assert that of the default.
    """
    from plasgenomicsutils.lib import ibd_selection as sel

    mat = _ibd_matrix(hotspot=900, seed=3)
    af = _af()
    got = {}
    for name, dt in (("f32", np.float32), ("f64", np.float64)):
        orig, sel._DTYPE = sel._DTYPE, dt
        try:
            got[name] = sel.compute_selection_statistic(mat, af, n_bins=10)[0]["raw_stat"]
        finally:
            sel._DTYPE = orig
    a, b = got["f32"], got["f64"]
    ok = np.isfinite(a) & np.isfinite(b)
    # same to single-precision relative tolerance, not merely correlated
    assert np.allclose(a[ok], b[ok], rtol=1e-4)
    assert np.nanmax(np.abs(b)) > 1.0          # O(1)+ , not a rounding residue


def test_the_published_variant_reproduces_the_cancellation_it_is_kept_for():
    from plasgenomicsutils.lib.ibd_selection import compute_selection_statistic

    mat = _ibd_matrix(hotspot=900, seed=3)
    af = _af()
    pub = compute_selection_statistic(mat, af, n_bins=10, variant="published")[0]
    cor = compute_selection_statistic(mat, af, n_bins=10, variant="corrected")[0]
    # The published recipe's per-SNP sum cancels, so all that is left is rounding residue --
    # tiny, and orders of magnitude below the corrected statistic. It keeps float32 (see
    # `_VARIANT_DTYPE`), so the residue sits near 1e-5 rather than 1e-14.
    assert np.nanmax(np.abs(pub["raw_stat"])) < 1e-3
    assert np.nanmax(np.abs(cor["raw_stat"])) > 1.0
    assert np.nanmax(np.abs(cor["raw_stat"])) / np.nanmax(np.abs(pub["raw_stat"])) > 1e4
    with pytest.raises(ValueError, match="variant must be"):
        compute_selection_statistic(mat, af, n_bins=10, variant="nope")


def test_the_upper_tail_stops_a_sharing_deficit_scoring_as_selection():
    from plasgenomicsutils.lib.ibd_selection import compute_selection_statistic

    # scattered segments, so the z-scores are symmetric noise and deficits exist at all: on
    # a real hotspot the corrected statistic runs one way only (about -0.2 to +10 here),
    # which is the point of it
    mat = _ibd_matrix(hotspot=None, seed=3)
    af = _af()
    two = compute_selection_statistic(mat, af, n_bins=10, tail="two-sided")[0]
    up = compute_selection_statistic(mat, af, n_bins=10, tail="upper")[0]
    assert np.allclose(two["z_score"], up["z_score"], equal_nan=True)   # same statistic

    deficit = np.isfinite(up["z_score"]) & (up["z_score"] < -2)
    assert deficit.any()
    # two-sided scores a deficit like an excess; the upper tail sends it to p -> 1
    assert np.nanmax(two["neg_log10_p"][deficit]) > 1
    assert np.nanmax(up["neg_log10_p"][deficit]) < 0.05
    # and the direction column names it either way
    assert set(np.asarray(up["direction"])[deficit]) == {"deficit"}
    excess = np.isfinite(up["z_score"]) & (up["z_score"] > 2)
    assert set(np.asarray(up["direction"])[excess]) == {"excess"}
    with pytest.raises(ValueError, match="tail must be"):
        compute_selection_statistic(mat, af, n_bins=10, tail="nope")
