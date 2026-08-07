"""Linkage-disequilibrium decay: mean r-squared against the distance between SNP pairs.

How far LD reaches is a property of the population -- in a high-transmission, freely
recombining setting it collapses within a few kb, and in a clonal or bottlenecked one it
persists. Read per metadata group, the curves separate populations that a single
genome-wide number would not.

This lives on the Python side because it is the one genuinely quadratic statistic in the
set: every SNP has to be compared with every other SNP within `max_dist`. The pairwise r
comes from :func:`allel.rogers_huff_r`, which counts each pair over the samples called at
both SNPs (pairwise-complete), so one missing call costs that pair rather than the SNP.

Two things bound the work. Each chromosome is thinned to `max_snps` evenly-spaced SNPs,
and only partners within `max_dist` are compared -- so the cost is linear in SNPs and
quadratic only in how many fall inside one `max_dist` window. Both are recorded in the
output header, because thinning that is invisible is thinning that gets forgotten.

Two conventions differ from the R implementation this replaced, both deliberately.
Distance bins are **half-open** `[start, end)`, matching the coordinate rule used
everywhere else in these packages (R's `cut()` closes on the right, which put a pair at
exactly a bin edge in the bin below). And thinning happens **after** the MAF filter, so
`max_snps` is a budget of usable SNPs rather than of candidates. On one chromosome with
thinning off, the two agree on all 19,459 pairs and on mean r-squared to 1e-5.

**r-squared is upward-biased in small samples**, by roughly 1/n. A group of ten sits above
a group of a hundred at every distance, and its half-decay can come back empty simply
because the curve never falls to half of an inflated first bin. Compare groups of similar
size, or subsample the larger ones.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: Pairs further apart than this are not counted.
LD_MAX_DIST = 50_000
#: SNPs kept per chromosome before the pairwise scan.
LD_MAX_SNPS = 3_000
#: Minor-allele-frequency floor; rare alleles give noisy, systematically low r-squared.
LD_MIN_MAF = 0.05
#: SNPs per block handed to rogers_huff_r at once.
_BLOCK = 2_000


def read_dosages(vcf_path, samples=None, regions=None, het="missing", min_depth=0):
    """Genotypes as a ``(n_variants, n_samples)`` int8 array for scikit-allel.

    Coded 0 = hom ref, 2 = hom alt, -1 = missing, which is what
    :func:`allel.rogers_huff_r` expects. The parasite is haploid, so ``het="missing"``
    (the default) reads a heterozygous call as a mixed infection and drops it at that
    site; ``het="dosage"`` keeps it as 1. Correlation is invariant to the 0/2 scaling, so
    this matches a 0/1 haploid coding exactly.

    Returns ``(gn, chrom, pos, samples)`` with ``pos`` 0-based.
    """
    from cyvcf2 import VCF

    vcf = VCF(vcf_path, samples=list(samples)) if samples else VCF(vcf_path)
    names = list(vcf.samples)
    if not names:
        raise SystemExit(f"{vcf_path}: no samples selected")

    rows, chroms, positions = [], [], []
    it = (v for r in regions for v in vcf(r)) if regions else vcf
    for v in it:
        gt = v.gt_types                       # 0 hom-ref, 1 het, 2 unknown, 3 hom-alt
        d = np.full(len(names), -1, dtype=np.int8)
        d[gt == 0] = 0
        d[gt == 3] = 2
        if het == "dosage":
            d[gt == 1] = 1
        if min_depth > 0:
            from .singletons import _read_depth
            depth = _read_depth(v, len(names))
            d[(depth >= 0) & (depth < min_depth)] = -1
        rows.append(d)
        chroms.append(v.CHROM)
        positions.append(v.POS - 1)           # VCF is 1-based; everything here is not
    vcf.close()
    if not rows:
        raise SystemExit(f"{vcf_path}: no variants read")
    return (np.vstack(rows), np.array(chroms), np.array(positions, dtype=np.int64), names)


def _maf(gn):
    """Minor-allele frequency per variant, ignoring missing calls."""
    called = gn >= 0
    n = called.sum(axis=1)
    alt = np.where(called, gn, 0).sum(axis=1) / 2.0     # 0/2 coding -> allele count
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(n > 0, alt / n, np.nan)
    return np.minimum(p, 1 - p), n


def _thin(idx, max_snps):
    """Keep at most `max_snps`, evenly spaced along the chromosome."""
    if len(idx) <= max_snps:
        return idx
    return idx[np.round(np.linspace(0, len(idx) - 1, max_snps)).astype(int)]


def _pairs_one_chrom(gn, pos, max_dist):
    """Yield ``(distance, r2)`` for every SNP pair on one chromosome within `max_dist`.

    rogers_huff_r returns the condensed matrix for *all* pairs it is given, so the
    chromosome is walked in blocks big enough to contain every pair inside `max_dist`,
    with each pair attributed to exactly one block.
    """
    import allel

    n = len(pos)
    if n < 2:
        return np.empty(0), np.empty(0)
    # furthest partner index still within max_dist, per SNP
    hi = np.searchsorted(pos, pos + max_dist, side="right") - 1
    reach = int(np.max(hi - np.arange(n)))
    if reach < 1:
        return np.empty(0), np.empty(0)

    d_out, r_out = [], []
    for start in range(0, n, _BLOCK):
        stop = min(n, start + _BLOCK + reach)
        sub = gn[start:stop]
        m = stop - start
        if m < 2:
            continue
        r = allel.rogers_huff_r(sub)                 # condensed, length m*(m-1)/2
        # condensed index of (i, j), i < j, within the block
        own = min(_BLOCK, m)                         # pairs owned by this block start at i < own
        for i in range(own):
            j_hi = min(hi[start + i] - start, m - 1)
            if j_hi <= i:
                continue
            j = np.arange(i + 1, j_hi + 1)
            k = (m * i - i * (i + 1) // 2) + (j - i - 1)
            d_out.append(pos[start + j] - pos[start + i])
            r_out.append(r[k])
        # no early exit: a block can reach the end of the chromosome while still owning
        # only its first `_BLOCK` SNPs, so the remaining ones need their own block
    if not d_out:
        return np.empty(0), np.empty(0)
    d = np.concatenate(d_out)
    r = np.concatenate(r_out)
    ok = np.isfinite(r) & (d <= max_dist)
    return d[ok], r[ok] ** 2


def ld_decay(gn, chrom, pos, groups=None, max_dist=LD_MAX_DIST, bins=20,
             maf=LD_MIN_MAF, max_snps=LD_MAX_SNPS, min_samples=4):
    """Binned r-squared against SNP-pair distance, per group.

    Args:
        gn: ``(n_variants, n_samples)`` int8 from :func:`read_dosages`.
        chrom, pos: per-variant chromosome and 0-based position.
        groups: per-sample group labels, or ``None`` to pool every sample.
        max_dist: longest pair separation to count, in bp.
        bins: number of distance bins, or an explicit sequence of edges.
        maf: minor-allele-frequency floor, applied within each group.
        max_snps: SNPs kept per chromosome before the scan.
        min_samples: skip groups smaller than this.

    Returns:
        ``(df, half)``: the binned table (`group`, `bin_start`, `bin_end`, `bin_mid`,
        `n_pairs`, `mean_r2`, `median_r2`) and a per-group half-decay table.
    """
    groups = np.array(["all"] * gn.shape[1]) if groups is None else np.asarray(groups)
    edges = (np.asarray(bins, dtype=float) if np.ndim(bins)
             else np.linspace(0, max_dist, int(bins) + 1))
    levels = [g for g in pd.unique(groups) if pd.notna(g)]

    rows = []
    for g in levels:
        cols = np.flatnonzero(groups == g)
        if len(cols) < min_samples:
            continue
        sub = gn[:, cols]
        m, _ = _maf(sub)
        dists, r2s = [], []
        for c in pd.unique(chrom):
            idx = np.flatnonzero((chrom == c) & (m >= maf) & np.isfinite(m))
            idx = _thin(idx, max_snps)
            if len(idx) < 2:
                continue
            d, r2 = _pairs_one_chrom(sub[idx], pos[idx], max_dist)
            if len(d):
                dists.append(d)
                r2s.append(r2)
        if not dists:
            continue
        d = np.concatenate(dists)
        r2 = np.concatenate(r2s)
        which = np.digitize(d, edges) - 1
        for b in range(len(edges) - 1):
            sel = which == b
            rows.append({
                "group": g, "bin_start": edges[b], "bin_end": edges[b + 1],
                "bin_mid": (edges[b] + edges[b + 1]) / 2,
                "n_pairs": int(sel.sum()),
                "mean_r2": float(r2[sel].mean()) if sel.any() else np.nan,
                "median_r2": float(np.median(r2[sel])) if sel.any() else np.nan,
            })
    df = pd.DataFrame(rows)
    return df, half_decay(df)


def half_decay(df):
    """Distance at which each group's binned mean first falls below half its first bin."""
    out = []
    for g, s in df.groupby("group", sort=False):
        s = s.dropna(subset=["mean_r2"]).sort_values("bin_mid")
        val = np.nan
        if len(s) >= 2:
            target = s["mean_r2"].iloc[0] / 2
            below = np.flatnonzero(s["mean_r2"].to_numpy() <= target)
            if len(below):
                k = int(below[0])
                if k == 0:
                    val = float(s["bin_mid"].iloc[0])
                else:
                    x0, x1 = s["bin_mid"].iloc[k - 1], s["bin_mid"].iloc[k]
                    y0, y1 = s["mean_r2"].iloc[k - 1], s["mean_r2"].iloc[k]
                    val = float(x1 if y0 == y1 else x0 + (target - y0) * (x1 - x0) / (y1 - y0))
        out.append({"group": g, "half_decay_bp": val})
    return pd.DataFrame(out)
