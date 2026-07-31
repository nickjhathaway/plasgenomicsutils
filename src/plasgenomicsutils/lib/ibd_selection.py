"""EIGENSTRAT-style IBD selection statistic (XiR,s).

Allele frequencies must be supplied externally and are never proxied from the
IBD matrix.

Method (Henden / Nygaard-style normalisation):
  1. binary IBD matrix (pairs x SNPs)
  2. subtract per-pair means, then per-SNP means
  3. divide by sqrt(p(1-p)), p = SNP allele frequency
  4. row-sum / sqrt(n_pairs) -> raw per-SNP statistic
  5. bin SNPs into equal-frequency MAF bins; within-bin z-score
  6. z^2 -> chi2(1df) -> p -> -log10(p); Bonferroni threshold
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2


# ---------------------------------------------------------------------------
# label / AF loading
# ---------------------------------------------------------------------------

def parse_snp_labels(snp_labels: list) -> pd.DataFrame:
    rows = []
    for label in snp_labels:
        if ":" in label:
            chrom, pos = label.rsplit(":", 1)
            rows.append({"snp_id": label, "chr": chrom, "pos": int(pos)})
        else:
            rows.append({"snp_id": label, "chr": "unknown", "pos": -1})
    return pd.DataFrame(rows)


def load_global_af(af_path: str, snp_labels: list) -> np.ndarray:
    """Global AFs aligned to ``snp_labels``; any missing SNP is a hard error."""
    af_df = pd.read_csv(af_path, sep="\t", usecols=["snp_id", "af"])
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
        print("  Hint: rerun allele_freqs with --zero-based to match 0-based matrix labels.")
        raise SystemExit(1)
    return af


def load_region_af_table(af_region_path: str) -> pd.DataFrame:
    df = pd.read_csv(af_region_path, sep="\t", usecols=["region", "snp_id", "af"])
    print(f"  Loaded region AF table: {len(df):,} rows, "
          f"{df['region'].nunique()} regions, {df['snp_id'].nunique():,} SNPs")
    return df


def get_af_for_region(region, snp_labels, region_af_table, global_af) -> np.ndarray:
    """Region AF with fallback: per-region table -> global AF. Never from the matrix."""
    if region_af_table is not None and region in region_af_table["region"].values:
        sub = region_af_table[region_af_table["region"] == region]
        af_map = sub.set_index("snp_id")["af"].to_dict()
        af = np.array([af_map.get(s, np.nan) for s in snp_labels])
        missing = int(np.isnan(af).sum())
        if missing > 0:
            print(f"  ERROR: {missing:,} SNPs (out of {len(snp_labels):,}) missing from "
                  f"region AF file for region '{region}'.")
            print("  Hint: rerun allele_freqs with --zero-based to match 0-based matrix labels.")
            raise SystemExit(1)
        return af
    return global_af.copy()


# ---------------------------------------------------------------------------
# pair subsetting
# ---------------------------------------------------------------------------

def within_region_row_indices(pair_labels, meta, region_col, region) -> np.ndarray:
    """Row indices for pairs where BOTH samples are from ``region``."""
    sample_to_region = meta.set_index("sample")[region_col].to_dict()
    idx = [
        i for i, label in enumerate(pair_labels)
        if all(sample_to_region.get(s) == region for s in label.split("__", 1))
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


def assemble_output(snp_df, stats, alpha, region=None) -> tuple[pd.DataFrame, float]:
    n_valid = int(np.sum(~np.isnan(stats["chi2_stat"])))
    threshold = -np.log10(alpha / n_valid) if n_valid > 0 else np.nan
    out = snp_df.copy()
    if region is not None:
        out.insert(0, "region", region)
    out["maf"] = stats["maf"]
    out["bin_id"] = stats["bin_id"]
    out["raw_stat"] = stats["raw_stat"]
    out["z_score"] = stats["z_score"]
    out["chi2_stat"] = stats["chi2_stat"]
    out["pval"] = stats["pval"]
    out["neg_log10_p"] = stats["neg_log10_p"]
    out["significant"] = out["neg_log10_p"] >= threshold
    return out, threshold
