"""Multiple-testing corrections for the IBD selection statistic.

Both a family-wise (Bonferroni) and a false-discovery-rate (Benjamini-Hochberg) cutoff are
reported, along with the genomic inflation factor that says whether either can be believed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
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
