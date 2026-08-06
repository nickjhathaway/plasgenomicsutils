"""EIGENSTRAT-style IBD selection statistic (XiR,s).

Allele frequencies must be supplied externally and are never proxied from the
IBD matrix.

Method (Henden / Nygaard-style normalisation):
  1. binary IBD matrix (pairs x SNPs)
  2. subtract per-pair means, then per-SNP means
  3. divide by sqrt(p(1-p)), p = SNP allele frequency
  4. row-sum / sqrt(n_pairs) -> raw per-SNP statistic
  5. bin SNPs into equal-frequency MAF bins; within-bin z-score
  6. z^2 -> chi2(1df) -> p -> -log10(p)
  7. two multiple-testing views: Bonferroni (family-wise) and Benjamini-Hochberg (FDR),
     plus the genomic inflation factor that says how much either can be trusted
  8. optionally `permutation_null`, which replaces steps 6-7's assumed reference with one
     drawn from the data

**On the thresholds.** Bonferroni asks for near-certainty that no SNP called is a false
positive; BH allows a stated share of them and so calls more. Both are reported, because
neither is obviously right here and the difference is worth seeing.

What is worth more than the choice between them is `lambda_gc`. The z-scores are
standardised to zero mean and unit variance *within each MAF bin*, which fixes the first
two moments and leaves the shape alone -- and the shape of IBD sharing is nothing like a
normal. On real *P. falciparum* data lambda comes out near 0.1 rather than 1: a tight bulk
with very heavy tails. Every p-value here descends from a chi2(1) that does not fit, so
both thresholds inherit that.

Rescaling cannot repair it. Standardising by median and MAD instead of mean and sd drives
lambda to exactly 1 by construction -- the diagnostic stops diagnosing -- while the tail
that actually decides significance gets further from normal, not closer.

`permutation_null` sidesteps the assumption. It re-draws the null by sliding each pair's
IBD segments to a random circular offset, so per-pair sharing, segment lengths and
along-genome autocorrelation all survive and only the alignment of segments *between*
pairs is destroyed. The 95th percentile of the per-replicate genome-wide maxima is then a
family-wise threshold that owes nothing to chi2(1), and the full null gives per-SNP
empirical p-values that a Benjamini-Hochberg pass can legitimately be run over. It costs
one full scan per replicate.

Note what that does and does not repair. `pval` and `q_value` still come from the chi2(1)
and remain miscalibrated -- do not quote them as probabilities. `p_empirical` and
`q_empirical` are the calibrated pair. The *ranking* was never in doubt either way: every
step from `z_score` to `neg_log10_p` is monotone, so the misfit moves no SNP relative to
another and only mislabels the axis.

Note also that no correction here knows about linkage. Adjacent SNPs in one sweep are not
independent tests, so the SNP counts overstate the number of findings; merge significant
SNPs into peaks before counting discoveries.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .intervals import SNP_COORD_SYSTEM, check_snp_coord_system
from ..utils.small_utils import Utils


# ---------------------------------------------------------------------------
# label / AF loading
# ---------------------------------------------------------------------------

def parse_snp_labels(snp_labels: list, with_pos_vcf: bool = False) -> pd.DataFrame:
    """Split ``chr:pos0`` labels into columns. ``pos`` is 0-based, like the label.

    ``with_pos_vcf`` adds the 1-based VCF position as ``pos_vcf`` for cross-referencing a
    variant against the VCF or a browser; it is redundant with ``pos``, so it is off by
    default rather than inflating every per-SNP table.
    """
    rows = []
    for label in snp_labels:
        if ":" in label:
            chrom, pos = label.rsplit(":", 1)
            rows.append({"snp_id": label, "chr": chrom, "pos": int(pos)})
        else:
            rows.append({"snp_id": label, "chr": "unknown", "pos": -1})
    df = pd.DataFrame(rows)
    if with_pos_vcf and len(df):
        df["pos_vcf"] = np.where(df["pos"] >= 0, df["pos"] + 1, -1)
    return df


def _read_af_table(af_path: str, usecols: list) -> pd.DataFrame:
    """Read an AF table, verifying its stamped SNP coordinate system first."""
    with Utils.smart_open_read(af_path) as fh:
        first = fh.readline().strip()
    stamped = first.split("=", 1)[1] if first.startswith("#snp_coord_system=") else None
    check_snp_coord_system(stamped, af_path)
    return pd.read_csv(af_path, sep="\t", usecols=usecols, comment="#")


def load_global_af(af_path: str, snp_labels: list) -> np.ndarray:
    """Global AFs aligned to ``snp_labels``; any missing SNP is a hard error."""
    af_df = _read_af_table(af_path, ["snp_id", "af"])
    af_map = af_df.set_index("snp_id")["af"].to_dict()
    af = np.array([af_map.get(s, np.nan) for s in snp_labels])
    missing = int(np.isnan(af).sum())
    if missing > 0:
        print(f"  ERROR: {missing:,} SNPs are missing from the AF file.")
        print("  First 10 SNP IDs expected (from matrix label file):")
        for s in snp_labels[:10]:
            print(f"    '{s}'")
        print("  First 10 SNP IDs found in AF file:")
        for s in list(af_map.keys())[:10]:
            print(f"    '{s}'")
        raise SystemExit(1)
    return af


def load_group_af_table(af_group_path: str) -> pd.DataFrame:
    df = _read_af_table(af_group_path, ["group", "snp_id", "af"])
    print(f"  Loaded group AF table: {len(df):,} rows, "
          f"{df['group'].nunique()} groups, {df['snp_id'].nunique():,} SNPs")
    return df


def get_af_for_group(group, snp_labels, group_af_table, global_af) -> np.ndarray:
    """Group AF with fallback: per-group table -> global AF. Never from the matrix."""
    if group_af_table is not None and group in group_af_table["group"].values:
        sub = group_af_table[group_af_table["group"] == group]
        af_map = sub.set_index("snp_id")["af"].to_dict()
        af = np.array([af_map.get(s, np.nan) for s in snp_labels])
        missing = int(np.isnan(af).sum())
        if missing > 0:
            print(f"  ERROR: {missing:,} SNPs (out of {len(snp_labels):,}) missing from "
                  f"group AF file for group '{group}'.")
            raise SystemExit(1)
        return af
    return global_af.copy()


# ---------------------------------------------------------------------------
# pair subsetting
# ---------------------------------------------------------------------------

def within_group_row_indices(pair_labels, meta, group_col, group) -> np.ndarray:
    """Row indices for pairs where BOTH samples are from ``group``."""
    sample_to_group = meta.set_index("sample")[group_col].to_dict()
    idx = [
        i for i, label in enumerate(pair_labels)
        if all(sample_to_group.get(s) == group for s in label.split("__", 1))
    ]
    return np.array(idx, dtype=np.int64)


# ---------------------------------------------------------------------------
# core statistic
# ---------------------------------------------------------------------------

def compute_selection_statistic(mat, af: np.ndarray, n_bins: int = 100,
                                label: str = "") -> tuple[dict, pd.DataFrame]:
    n_pairs, n_snps = mat.shape
    tag = f"[{label}] " if label else ""
    print(f"  {tag}n_pairs={n_pairs:,}  n_snps={n_snps:,}")
    estimated_gb = n_pairs * n_snps * 4 / 1e9
    if estimated_gb > 16:
        print(f"  {tag}~{estimated_gb:.1f} GB dense — using chunked path")
        return _compute_chunked(mat, af, n_bins)
    print(f"  {tag}~{estimated_gb:.1f} GB dense — using dense path")
    return _compute_dense(mat, af, n_bins)


def _compute_dense(mat, af, n_bins):
    n_pairs, n_snps = mat.shape
    X = mat.toarray().astype(np.float32).T  # (snps, pairs)
    X -= X.mean(axis=0, keepdims=True)      # per-pair means
    X -= X.mean(axis=1, keepdims=True)      # per-SNP means
    valid = ~np.isnan(af) & (af > 0) & (af < 1)
    denom = np.where(valid, np.sqrt(af * (1 - af)), np.nan).astype(np.float32)
    X /= denom[:, np.newaxis]
    raw_stat = np.nansum(X, axis=1) / np.sqrt(n_pairs)
    raw_stat = np.where(valid, raw_stat, np.nan)
    return _normalise_and_finalise(raw_stat, af, valid, n_bins)


def _compute_chunked(mat, af, n_bins, chunk_size=500):
    n_pairs, n_snps = mat.shape
    pair_means = np.asarray(mat.mean(axis=1)).ravel().astype(np.float32)
    valid = ~np.isnan(af) & (af > 0) & (af < 1)
    raw_stat = np.full(n_snps, np.nan, dtype=np.float64)
    n_chunks = (n_snps + chunk_size - 1) // chunk_size
    for c in range(n_chunks):
        start = c * chunk_size
        end = min(start + chunk_size, n_snps)
        if c % max(1, n_chunks // 10) == 0:
            print(f"    chunk {c+1}/{n_chunks}  SNPs {start}-{end}", end="\r", flush=True)
        chunk = mat[:, start:end].toarray().astype(np.float32).T  # (chunk, pairs)
        chunk -= pair_means[np.newaxis, :]
        chunk -= chunk.mean(axis=1, keepdims=True)
        chunk_af = af[start:end]
        chunk_valid = valid[start:end]
        denom = np.where(chunk_valid, np.sqrt(chunk_af * (1 - chunk_af)), np.nan).astype(np.float32)
        chunk /= denom[:, np.newaxis]
        rs = np.nansum(chunk, axis=1) / np.sqrt(n_pairs)
        raw_stat[start:end] = np.where(chunk_valid, rs, np.nan)
    print()
    return _normalise_and_finalise(raw_stat, af, valid, n_bins)


def _normalise_and_finalise(raw_stat, af, valid, n_bins):
    n_snps = len(raw_stat)
    maf = np.where(af <= 0.5, af, 1 - af)

    bin_ids = np.full(n_snps, -1, dtype=int)
    valid_idx = np.where(valid)[0]
    maf_valid = maf[valid_idx]
    bin_edges = np.quantile(maf_valid, np.linspace(0, 1, n_bins + 1))
    bin_edges[-1] += 1e-9
    bin_ids[valid_idx] = np.digitize(maf_valid, bin_edges) - 1

    z_score = np.full(n_snps, np.nan)
    bin_records = []
    for b in range(n_bins):
        idx = np.where(bin_ids == b)[0]
        if len(idx) < 2:
            continue
        vals = raw_stat[idx]
        mu = np.nanmean(vals)
        sd = np.nanstd(vals, ddof=1)
        if sd == 0 or np.isnan(sd):
            continue
        z_score[idx] = (vals - mu) / sd
        bin_records.append({
            "bin": b, "n_snps": len(idx),
            "maf_min": bin_edges[b], "maf_max": bin_edges[b + 1],
            "mean": mu, "sd": sd,
        })
    bin_df = pd.DataFrame(bin_records)

    chi2_stat = z_score ** 2
    pval = np.where(~np.isnan(chi2_stat), chi2.sf(chi2_stat, df=1), np.nan)
    neg_log10_p = np.where(pval > 0, -np.log10(pval), np.nan)

    return {
        "raw_stat": raw_stat, "z_score": z_score,
        "chi2_stat": chi2_stat, "pval": pval,
        "neg_log10_p": neg_log10_p, "bin_id": bin_ids, "maf": maf,
    }, bin_df


def benjamini_hochberg(pval):
    """Benjamini-Hochberg q-values; `NaN` in, `NaN` out.

    BH controls the expected *proportion* of false positives among the SNPs called, where
    Bonferroni controls the probability of even one. It is valid under independence and
    under positive regression dependency -- the usual justification for using it on a
    genome scan, since SNPs in linkage are positively correlated. Treat that as an
    approximation: BH is now known not to control FDR in general for correlated two-sided
    tests. The larger practical worry is upstream of either correction, in whether the
    p-values are calibrated at all -- see :func:`genomic_inflation`.
    """
    p = np.asarray(pval, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.flatnonzero(~np.isnan(p))
    if not ok.size:
        return q
    m = ok.size
    order = ok[np.argsort(p[ok])]
    ranked = p[order] * m / np.arange(1, m + 1)
    q[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    return q


def genomic_inflation(chi2_stat):
    """Median chi2 over its null expectation -- 1.0 when the reference is right.

    Far from 1 means the chi2(1) null is the wrong distribution, and every p-value drawn
    from it is wrong with it. That happens easily here: the z-scores are standardised to
    unit variance within each MAF bin, which pins the first two moments and says nothing
    about the shape, while IBD sharing is autocorrelated and heavy-tailed. Read this
    before reading any threshold.
    """
    v = np.asarray(chi2_stat, dtype=float)
    v = v[~np.isnan(v)]
    if not v.size:
        return np.nan
    return float(np.median(v) / chi2.ppf(0.5, 1))


def permutation_null(mat, af, n_perm=200, n_bins=100, alpha=0.05, seed=0,
                     progress=None):
    """Null distribution for the selection statistic, built by moving the IBD around.

    Neither Bonferroni nor Benjamini-Hochberg is trustworthy on this statistic. Both
    assume the chi2(1) reference fits, and it does not (see :func:`genomic_inflation`),
    and both treat SNPs as independent tests when one IBD segment spans hundreds of them.
    This generates the null instead of assuming it.

    Each replicate slides **every pair's IBD segments to a random position** along the SNP
    axis, wrapping at the end. That keeps each pair's total sharing, its segment count and
    its segment lengths exactly as observed -- so relatedness, block structure and the
    resulting autocorrelation all survive -- and destroys only the thing being tested:
    whether pairs share *the same* locus. The statistic is then recomputed from scratch,
    MAF binning and all.

    Three summaries come out of the one pass, in increasing order of resolution and of
    what they assume:

    ``threshold``
        The ``1 - alpha`` quantile of the per-replicate genome-wide **maxima**: the
        largest score reachable with no locus-specific sharing. A family-wise line, and
        the assumption-free one.
    ``p_pointwise``
        Per SNP, how often the null at *that same SNP* reached the observed value. Assumes
        nothing beyond the shift itself, but cannot resolve below ``1 / (n_perm + 1)``, so
        it is a sanity check on the top hits rather than an input to FDR.
    ``p_pooled``
        Per SNP, how often *any* null value anywhere in the genome reached the observed
        value. Resolution ``1 / (n_perm * n_snps + 1)`` -- fine enough to feed
        Benjamini-Hochberg -- at the price of assuming the null is exchangeable across
        SNPs. The within-MAF-bin standardisation is what is supposed to make that true;
        ``bin_tail_rate`` reports whether it did (see below).
    ``p_stratified``
        The same count taken only within the SNP's own MAF bin, which drops the
        exchangeability assumption entirely. Correct where ``p_pooled`` is approximate,
        but ``n_bins`` times coarser -- with the defaults that is a resolution of about
        1e-5, too blunt for a genome-wide FDR, so it is a cross-check rather than a
        replacement unless ``n_perm`` is raised by roughly the same factor.

    Both p-values use the Phipson-Smyth ``(1 + exceedances) / (1 + draws)`` form, so no
    SNP is ever assigned p = 0 -- a permutation cannot evidence a p smaller than its own
    resolution.

    ``bin_tail_rate`` is the exchangeability check behind ``p_pooled``: the share of each
    MAF bin's null values landing above one common reference (the 99th percentile of the
    first replicate). Under exchangeability every bin returns 0.01; a bin far from that is
    a bin whose null is shaped differently from the rest, and pooling across it is doing
    the SNPs in it a quiet disservice.

    Expect the family-wise line to be far stricter than the parametric ones. Measured on
    one group of a real 249-sample cohort (4,278 pairs, 27,897 SNPs, 200 replicates), the
    null maxima ran 7.4 to 33.8 with a median of 12.7, against a Bonferroni line of 5.73
    -- so **every one of the 200 replicates produced a genome-wide maximum that Bonferroni
    would have called significant**, on data constructed to hold no shared locus at all. A
    nominal 5% family-wise line with an actual family-wise error rate of 100% is the
    concrete cost of the chi2(1) misfit.

    Args:
        mat: The (pairs x SNPs) binary matrix for one group.
        af: Allele frequencies aligned to its columns.
        n_perm: Replicates. 200 is enough for a 5% quantile and gives ``p_pooled`` about
            six significant figures of headroom; raise it for a smaller alpha.
        n_bins, alpha, seed: MAF bins, family-wise level, RNG seed.
        progress: Optional ``callable(i, n_perm)`` for a progress line.

    Returns:
        A dict with ``threshold``, ``maxima``, ``p_pointwise``, ``p_pooled``,
        ``bin_tail_rate``, ``n_perm`` and ``n_pool``. Everything per-SNP is aligned to the
        columns of ``mat`` and carries ``NaN`` wherever the observed statistic does.
    """
    from scipy import sparse

    coo = mat.tocoo()
    n_pairs, n_snps = mat.shape
    obs = compute_selection_statistic(mat, af, n_bins=n_bins, label="")[0]["neg_log10_p"]
    obs_ok = np.isfinite(obs)

    rng = np.random.default_rng(seed)
    maxima = np.empty(n_perm)
    point = np.zeros(n_snps, dtype=np.int64)     # null beat the observed AT this SNP
    pooled = np.zeros(n_snps, dtype=np.int64)    # null beat it ANYWHERE
    strat = np.zeros(n_snps, dtype=np.int64)     # ...anywhere IN ITS OWN MAF BIN
    bin_n = np.zeros(n_bins, dtype=np.int64)
    bin_hi = np.zeros(n_bins, dtype=np.int64)
    n_pool = 0
    ref = np.nan
    # Bin membership is a function of `af` and `n_bins` alone, so it is the same in every
    # replicate as in the observed scan -- which is what lets the stratified count below
    # compare a SNP only against nulls from its own bin.
    members = None

    for r in range(n_perm):
        shift = rng.integers(0, n_snps, size=n_pairs)
        col = (coo.col + shift[coo.row]) % n_snps
        null = sparse.coo_matrix((coo.data, (coo.row, col)), shape=mat.shape).tocsr()
        st, _ = compute_selection_statistic(null, af, n_bins=n_bins, label="")
        v = st["neg_log10_p"]
        maxima[r] = np.nanmax(v)

        # NaN compares False, so unusable SNPs simply never count as exceedances
        point += v >= obs
        good = np.sort(v[np.isfinite(v)])
        n_pool += good.size
        # searchsorted puts NaN past the end, giving those SNPs a count of 0; they are
        # masked to NaN below rather than silently reported as significant
        pooled += good.size - np.searchsorted(good, obs, side="left")

        b = st["bin_id"]
        keep = np.isfinite(v) & (b >= 0)
        if r == 0:
            ref = float(np.nanquantile(v, 0.99))
            members = [np.where(b == k)[0] for k in range(n_bins)]
        bin_n += np.bincount(b[keep], minlength=n_bins)[:n_bins]
        bin_hi += np.bincount(b[keep & (v >= ref)], minlength=n_bins)[:n_bins]

        for idx in members:
            if idx.size == 0:
                continue
            vb = np.sort(v[idx][np.isfinite(v[idx])])
            if vb.size:
                strat[idx] += vb.size - np.searchsorted(vb, obs[idx], side="left")

        if progress is not None:
            progress(r + 1, n_perm)

    # draws behind each SNP's stratified p: its own bin's finite nulls, over all replicates
    n_strat = np.zeros(n_snps, dtype=np.int64)
    for k, idx in enumerate(members or []):
        n_strat[idx] = bin_n[k]

    p_point = np.where(obs_ok, (1 + point) / (1 + n_perm), np.nan)
    p_pool = np.where(obs_ok, (1 + pooled) / (1 + n_pool), np.nan)
    p_strat = np.where(obs_ok & (n_strat > 0),
                       (1 + strat) / (1 + np.maximum(n_strat, 1)), np.nan)
    # Only bins with enough draws to resolve a 1% rate get one: below ~1000 the estimate
    # is dominated by counting noise, and a small bin returning 0 would otherwise make the
    # spread look infinite when nothing is actually wrong with it.
    _MIN_DRAWS = 1000
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = np.where(bin_n >= _MIN_DRAWS, bin_hi / np.maximum(bin_n, 1), np.nan)

    return {
        "threshold": float(np.quantile(maxima, 1 - alpha)),
        "maxima": maxima,
        "p_pointwise": p_point,
        "p_pooled": p_pool,
        "p_stratified": p_strat,
        "n_stratified": n_strat,
        "bin_tail_rate": rate,
        "bin_tail_rate_skipped": int(np.sum((bin_n > 0) & (bin_n < _MIN_DRAWS))),
        "n_perm": int(n_perm),
        "n_pool": int(n_pool),
    }


def assemble_output(snp_df, stats, alpha, group=None, fdr_alpha=None, perm=None,
                    pool="global"):
    """Per-SNP table plus every threshold and the calibration diagnostic.

    ``perm`` is the dict from :func:`permutation_null`, or ``None`` to skip the
    permutation columns entirely. ``pool`` chooses which of its empirical p-values drives
    ``q_empirical``: ``"global"`` pools the null across all SNPs (fine resolution, assumes
    exchangeability) and ``"bin"`` stays inside each MAF bin (assumption-free, ``n_bins``
    times coarser). Both are written either way; only the q-value follows ``pool``.

    Returns ``(df, info)``. ``info`` holds each cutoff on the ``-log10(p)`` scale so any
    of them can be drawn as a line, and ``lambda_gc``, which says how much the parametric
    ones are worth.
    """
    fdr_alpha = alpha if fdr_alpha is None else fdr_alpha
    n_valid = int(np.sum(~np.isnan(stats["chi2_stat"])))
    threshold = -np.log10(alpha / n_valid) if n_valid > 0 else np.nan

    out = snp_df.copy()
    if group is not None:
        out.insert(0, "group", group)
    out["maf"] = stats["maf"]
    out["bin_id"] = stats["bin_id"]
    out["raw_stat"] = stats["raw_stat"]
    out["z_score"] = stats["z_score"]
    out["chi2_stat"] = stats["chi2_stat"]
    out["pval"] = stats["pval"]
    out["neg_log10_p"] = stats["neg_log10_p"]
    out["q_value"] = benjamini_hochberg(stats["pval"])
    out["significant"] = out["neg_log10_p"] >= threshold          # Bonferroni (FWER)
    out["significant_fdr"] = out["q_value"] < fdr_alpha           # Benjamini-Hochberg
    perm_threshold = perm["threshold"] if perm else None
    if perm_threshold is not None and np.isfinite(perm_threshold):
        out["significant_perm"] = out["neg_log10_p"] >= perm_threshold
    if perm:
        # Calibrated companions to `pval`/`q_value`: same SNPs, same order, but a
        # reference drawn from the data rather than assumed. `p_pointwise` bottoms out at
        # 1/(n_perm+1) so it cannot drive FDR; it is here as the assumption-free check on
        # whatever the pooled p-value calls.
        if pool not in ("global", "bin"):
            raise ValueError("pool must be 'global' or 'bin'")
        out["p_pointwise"] = perm["p_pointwise"]
        out["p_empirical"] = perm["p_pooled"]
        out["p_empirical_binned"] = perm["p_stratified"]
        chosen = perm["p_pooled"] if pool == "global" else perm["p_stratified"]
        out["q_empirical"] = benjamini_hochberg(chosen)
        out["significant_fdr_perm"] = out["q_empirical"] < fdr_alpha

    # The BH critical value for the number of rejections, k * q / m, as -log10(p) so it
    # plots as a line. Not the largest p actually called: this is the cutoff BH applies,
    # so `neg_log10_p >= line` reproduces the flag exactly. It is always at or below the
    # Bonferroni line (which is the k = 1 case).
    n_rej = int(out["significant_fdr"].fillna(False).sum())
    fdr_threshold = (float(-np.log10(n_rej * fdr_alpha / n_valid))
                     if n_rej and n_valid else np.nan)

    # The smallest q the permutation could possibly produce: the top SNP's p bottoms out
    # at 1/(1 + n_perm * n_snps), and BH multiplies the rank-1 p-value by m. So the whole
    # empirical-FDR column is dead unless this sits below `fdr_alpha` -- and it wants to
    # sit well below, since one stray null exceedance at the top SNP multiplies it. In
    # round terms q_floor ~ 1 / n_perm, so n_perm >= 10 / fdr_alpha buys an order of
    # magnitude of headroom (200 replicates at the default q < 0.05).
    # The smallest q each SNP could reach: BH multiplies its p by the number of tests, and
    # its p bottoms out at 1/(1 + draws behind it). Under global pooling every SNP has the
    # same draws and so the same floor; under bin pooling the floor is per-SNP, and MAF
    # bins are wildly unequal in practice (ties in MAF collapse the quantile edges), so a
    # single "best case" number would hide most of the column being unreachable.
    # How the MAF binning actually came out. Requesting N equal-frequency bins does not
    # give N: allele frequency is k/n for a smallish n, so MAF ties collapse the quantile
    # edges. This is a property of the statistic, not of the permutation, so it is
    # computed either way -- the within-bin standardisation is step 5 of the method.
    usable = ~np.isnan(stats["chi2_stat"])
    counts = np.bincount(stats["bin_id"][usable & (stats["bin_id"] >= 0)])
    n_bins_used = int((counts > 0).sum())
    big_bin = float(counts.max() / n_valid) if counts.size and n_valid else np.nan

    q_floor, frac_dead = np.nan, np.nan
    if perm and n_valid:
        draws = (np.full(len(stats["bin_id"]), perm["n_pool"]) if pool == "global"
                 else perm["n_stratified"])
        with np.errstate(invalid="ignore", divide="ignore"):
            per_snp = np.where(draws > 0, n_valid / (1.0 + draws), np.inf)
        q_floor = float(np.min(per_snp[usable]))
        frac_dead = float(np.mean(per_snp[usable] > fdr_alpha))

    # The empirical-FDR line, as the smallest observed score BH kept. Unlike the k*q/m
    # form above there is no closed expression for it, because the empirical p-values are
    # a step function of the score rather than a smooth transform of it.
    emp_line = np.nan
    if perm and out["significant_fdr_perm"].any():
        emp_line = float(out.loc[out["significant_fdr_perm"].fillna(False),
                                 "neg_log10_p"].min())

    info = {
        "alpha": alpha, "n_tests": n_valid, "neg_log10_p_threshold": threshold,
        "neg_log10_p_perm_threshold": (float(perm_threshold)
                                       if perm_threshold is not None else np.nan),
        "n_significant_perm": (int(out["significant_perm"].sum())
                               if "significant_perm" in out else -1),
        "fdr_alpha": fdr_alpha, "neg_log10_p_fdr_threshold": fdr_threshold,
        "n_significant": int(out["significant"].sum()),
        "n_significant_fdr": int(out["significant_fdr"].fillna(False).sum()),
        "neg_log10_p_emp_fdr_threshold": emp_line,
        "n_significant_fdr_perm": (int(out["significant_fdr_perm"].fillna(False).sum())
                                   if perm else -1),
        "n_perm": perm["n_perm"] if perm else 0,
        "p_empirical_resolution": (1.0 / (1 + perm["n_pool"])) if perm else np.nan,
        "q_empirical_floor": q_floor, "empirical_pool": (pool if perm else ""),
        "frac_q_unreachable": frac_dead,
        "n_bins_used": n_bins_used, "largest_bin_frac": big_bin,
        "perm_bin_tail_min": (float(np.nanmin(perm["bin_tail_rate"])) if perm else np.nan),
        "perm_bin_tail_max": (float(np.nanmax(perm["bin_tail_rate"])) if perm else np.nan),
        "lambda_gc": genomic_inflation(stats["chi2_stat"]),
    }
    return out, info
