"""LD decay. The heavy pairwise work is scikit-allel's rogers_huff_r; these tests pin the
windowing, binning and grouping around it, and that missing calls cost a pair rather than
a SNP."""

from __future__ import annotations

import numpy as np
import pytest

from plasgenomicsutils.lib.ld import _pairs_one_chrom, _thin, half_decay, ld_decay


def _linked(n_samples=60, n_var=12, noise=0.0, seed=0):
    """Variants along a chromosome, each copying its neighbour with a per-step flip
    probability. Correlation then decays as (1 - 2*noise)^steps, so r-squared genuinely
    falls with distance -- a flip probability that grows past 0.5 instead produces
    *anti*-correlation and r-squared climbing back up, which is not decay."""
    rng = np.random.default_rng(seed)
    cur = rng.integers(0, 2, n_samples)
    rows = [cur.copy()]
    for _ in range(n_var - 1):
        cur = np.where(rng.random(n_samples) < noise, 1 - cur, cur)
        rows.append(cur.copy())
    return np.array(rows, dtype=np.int8) * 2


def test_perfectly_linked_variants_have_r2_of_one():
    gn = _linked(noise=0.0)
    pos = np.arange(gn.shape[0]) * 1000
    d, r2 = _pairs_one_chrom(gn, pos, max_dist=50_000)
    assert len(d) == gn.shape[0] * (gn.shape[0] - 1) // 2      # every pair is in range
    assert r2 == pytest.approx(1.0)


def test_only_pairs_within_max_dist_are_counted():
    gn = _linked(n_var=6)
    pos = np.array([0, 1000, 2000, 50_000, 51_000, 52_000])
    d, _ = _pairs_one_chrom(gn, pos, max_dist=5000)
    # three pairs inside each cluster, none across the 48 kb gap
    assert len(d) == 6
    assert d.max() <= 5000


def test_blocking_does_not_change_the_pair_set(monkeypatch):
    """Blocks exist only to bound memory; each pair must be counted exactly once."""
    import plasgenomicsutils.lib.ld as ld

    gn = _linked(n_var=40, noise=0.01)
    pos = np.arange(gn.shape[0]) * 1000
    d_ref, r_ref = _pairs_one_chrom(gn, pos, max_dist=10_000)
    monkeypatch.setattr(ld, "_BLOCK", 5)                        # force many small blocks
    d_small, r_small = ld._pairs_one_chrom(gn, pos, max_dist=10_000)
    order_ref, order_small = np.argsort(d_ref, kind="stable"), np.argsort(d_small, kind="stable")
    assert len(d_ref) == len(d_small)
    assert np.array_equal(d_ref[order_ref], d_small[order_small])
    assert r_ref[order_ref] == pytest.approx(r_small[order_small], abs=1e-6)


def test_a_missing_call_costs_the_pair_not_the_snp():
    gn = _linked(n_var=4)
    pos = np.arange(4) * 1000
    full = _pairs_one_chrom(gn, pos, max_dist=50_000)[1]
    gm = gn.copy()
    gm[0, :20] = -1                                             # a third of one SNP
    part = _pairs_one_chrom(gm, pos, max_dist=50_000)[1]
    assert len(part) == len(full)                               # the SNP is still scanned
    assert np.all(np.isfinite(part))


def test_decay_falls_with_distance_and_bins_are_half_open():
    gn = _linked(n_var=40, noise=0.05, seed=3)
    pos = np.arange(gn.shape[0]) * 1000
    chrom = np.array(["c1"] * gn.shape[0])
    df, half = ld_decay(gn, chrom, pos, max_dist=40_000, bins=4, maf=0.0)
    ok = df.dropna(subset=["mean_r2"])
    assert ok["mean_r2"].iloc[0] > ok["mean_r2"].iloc[-1]
    # bins tile [0, max_dist) with no gaps and no overlap
    assert list(ok["bin_start"]) == [0, 10_000, 20_000, 30_000]
    assert list(ok["bin_end"]) == [10_000, 20_000, 30_000, 40_000]
    assert int(df["n_pairs"].sum()) == len(_pairs_one_chrom(gn, pos, 40_000)[0])


def test_groups_are_scanned_separately_and_small_ones_skipped():
    gn = np.hstack([_linked(n_samples=40, n_var=10, noise=0.0),          # group a: linked
                    _linked(n_samples=40, n_var=10, noise=0.5, seed=7),  # group b: not
                    _linked(n_samples=2, n_var=10)])                     # too small
    pos = np.arange(10) * 1000
    chrom = np.array(["c1"] * 10)
    groups = np.array(["a"] * 40 + ["b"] * 40 + ["tiny"] * 2)
    df, _ = ld_decay(gn, chrom, pos, groups=groups, max_dist=20_000, bins=2, maf=0.0)
    assert set(df["group"]) == {"a", "b"}                       # 'tiny' is below min_samples
    a = df[df.group == "a"]["mean_r2"].iloc[0]
    b = df[df.group == "b"]["mean_r2"].iloc[0]
    assert a > b


def test_thinning_is_even_and_bounded():
    idx = np.arange(1000)
    assert np.array_equal(_thin(idx, 5000), idx)                # under budget: untouched
    t = _thin(idx, 100)
    assert len(t) == 100
    assert t[0] == 0 and t[-1] == 999                           # spans the chromosome
    assert np.all(np.diff(t) > 0)


def test_half_decay_interpolates_and_reports_nothing_when_flat():
    import pandas as pd

    df = pd.DataFrame({"group": ["a"] * 3 + ["b"] * 3,
                       "bin_mid": [1000, 2000, 3000] * 2,
                       "mean_r2": [0.4, 0.3, 0.1, 0.2, 0.2, 0.2]})
    h = half_decay(df).set_index("group")["half_decay_bp"]
    assert 2000 < h["a"] < 3000        # crosses 0.2 between the second and third bins
    assert np.isnan(h["b"])            # never halves
