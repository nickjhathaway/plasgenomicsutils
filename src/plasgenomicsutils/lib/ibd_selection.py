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

**On the thresholds.** Bonferroni asks for near-certainty that no SNP called is a false
positive; BH allows a stated share of them and so calls more. Both are reported, because
neither is obviously right here and the difference is worth seeing.

What is worth more than the choice between them is `lambda_gc`. The z-scores are
standardised to zero mean and unit variance *within each MAF bin*, which fixes the first
two moments and leaves the shape alone -- and the shape of IBD sharing is nothing like a
normal. On real *P. falciparum* data lambda comes out near 0.1 rather than 1: a tight bulk
with very heavy tails. Every p-value here descends from a chi2(1) that does not fit, so
both thresholds inherit that, and the honest use of either is as a *ranking* device.

Also note that neither correction knows about linkage. Adjacent SNPs in one sweep are not
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


def assemble_output(snp_df, stats, alpha, group=None, fdr_alpha=None):
    """Per-SNP table plus both thresholds and the calibration diagnostic.

    Returns ``(df, info)``. ``info`` holds the Bonferroni and Benjamini-Hochberg cutoffs
    on the ``-log10(p)`` scale, so either can be drawn as a line, and ``lambda_gc``.
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

    # The BH critical value for the number of rejections, k * q / m, as -log10(p) so it
    # plots as a line. Not the largest p actually called: this is the cutoff BH applies,
    # so `neg_log10_p >= line` reproduces the flag exactly. It is always at or below the
    # Bonferroni line (which is the k = 1 case).
    n_rej = int(out["significant_fdr"].fillna(False).sum())
    fdr_threshold = (float(-np.log10(n_rej * fdr_alpha / n_valid))
                     if n_rej and n_valid else np.nan)

    info = {
        "alpha": alpha, "n_tests": n_valid, "neg_log10_p_threshold": threshold,
        "fdr_alpha": fdr_alpha, "neg_log10_p_fdr_threshold": fdr_threshold,
        "n_significant": int(out["significant"].sum()),
        "n_significant_fdr": int(out["significant_fdr"].fillna(False).sum()),
        "lambda_gc": genomic_inflation(stats["chi2_stat"]),
    }
    return out, info
