"""The selection statistic's dense and chunked paths must agree.

`compute_selection_statistic` switches to a chunked implementation once the dense matrix
would exceed ~16 GB, so on any test-sized input only the dense branch runs. Calling the two
private paths directly is the only way to pin that they compute the same thing -- otherwise
a change to one silently diverges from the other on exactly the large cohorts that use it.

The assertions are on ``raw_stat`` (what the two implementations actually compute
differently) and on the MAF bin assignment. The downstream ``z_score`` / ``chi2_stat`` /
``neg_log10_p`` are a shared function of those two and are deliberately *not* compared
element-wise: the arithmetic is float32, so ``raw_stat`` agrees only to ~1e-7, and z-scoring
divides by a within-bin standard deviation that on signal-free synthetic data is itself
~1e-7 -- which inflates that harmless difference to order 1. On real data the within-bin
spread is far larger than float32 noise.
"""

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from plasgenomicsutils.lib.ibd_selection import (
    _compute_chunked,
    _compute_dense,
    compute_selection_statistic,
)

RAW_TOL = 1e-6          # float32 accumulation over the pair axis


def _fixture(n_pairs=60, n_snps=120, seed=0, sweep=False):
    rng = np.random.default_rng(seed)
    m = (rng.random((n_pairs, n_snps)) < 0.25).astype(np.uint8)
    if sweep:                       # a real IBD sweep, so raw_stat has genuine spread
        m[: int(n_pairs * 0.7), 40:60] = 1
    return csr_matrix(m), rng.uniform(0.05, 0.95, size=n_snps)


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("chunk_size", [1, 7, 17, 500])
def test_raw_stat_matches_the_dense_path_at_any_chunk_size(seed, chunk_size):
    mat, af = _fixture(seed=seed)
    dense = np.asarray(_compute_dense(mat, af, n_bins=10)[0]["raw_stat"], dtype=float)
    chunked = np.asarray(
        _compute_chunked(mat, af, n_bins=10, chunk_size=chunk_size)[0]["raw_stat"],
        dtype=float)
    np.testing.assert_allclose(chunked, dense, rtol=0, atol=RAW_TOL, equal_nan=True)


def test_bin_assignment_is_identical():
    mat, af = _fixture(n_pairs=200, n_snps=400, sweep=True)
    dense = _compute_dense(mat, af, n_bins=10)[0]
    chunked = _compute_chunked(mat, af, n_bins=10, chunk_size=17)[0]
    assert np.array_equal(np.asarray(dense["bin_id"]), np.asarray(chunked["bin_id"]))
    np.testing.assert_allclose(np.asarray(dense["maf"], dtype=float),
                               np.asarray(chunked["maf"], dtype=float))


def test_a_ragged_final_chunk_is_handled():
    # 91 SNPs in chunks of 7 divides evenly; in chunks of 8 the last chunk is short
    mat, af = _fixture(n_pairs=37, n_snps=91, seed=3)
    even = np.asarray(_compute_chunked(mat, af, n_bins=5, chunk_size=7)[0]["raw_stat"], float)
    ragged = np.asarray(_compute_chunked(mat, af, n_bins=5, chunk_size=8)[0]["raw_stat"], float)
    np.testing.assert_allclose(ragged, even, rtol=0, atol=RAW_TOL, equal_nan=True)


def test_unusable_allele_frequencies_are_excluded_on_both_paths():
    mat, af = _fixture(n_snps=40)
    af[:4] = [np.nan, 0.0, 1.0, -0.1]        # missing / no variance / out of range
    with np.errstate(invalid="ignore"):
        dense = _compute_dense(mat, af, n_bins=5)[0]
        chunked = _compute_chunked(mat, af, n_bins=5, chunk_size=6)[0]
    for stats in (dense, chunked):
        assert np.isnan(np.asarray(stats["raw_stat"], dtype=float)[:4]).all()
        assert (np.asarray(stats["bin_id"])[:4] == -1).all()


def test_small_input_takes_the_dense_path(capsys):
    mat, af = _fixture(n_pairs=10, n_snps=20)
    compute_selection_statistic(mat, af, n_bins=4)
    assert "dense path" in capsys.readouterr().out
