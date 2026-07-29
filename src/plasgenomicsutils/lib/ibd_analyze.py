"""Downstream analysis of the binary IBD matrix.

Per-pair / per-SNP / per-region / per-chromosome IBD summaries.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .ibd_selection import parse_snp_labels


def parse_pair_labels(pair_labels: list) -> pd.DataFrame:
    rows = []
    for label in pair_labels:
        s1, s2 = label.split("__", 1)
        rows.append({"pair": label, "sample1": s1, "sample2": s2})
    return pd.DataFrame(rows)


def per_pair_summary(mat, pair_labels: list) -> pd.DataFrame:
    n_snps = mat.shape[1]
    counts = np.asarray(mat.sum(axis=1)).ravel()
    pair_df = parse_pair_labels(pair_labels)
    pair_df["n_ibd_snps"] = counts
    pair_df["frac_ibd"] = counts / n_snps
    return pair_df


def per_snp_summary(mat, snp_labels: list) -> pd.DataFrame:
    n_pairs = mat.shape[0]
    counts = np.asarray(mat.sum(axis=0)).ravel()
    snp_df = parse_snp_labels(snp_labels)
    snp_df["n_pairs_ibd"] = counts
    snp_df["frac_pairs_ibd"] = counts / n_pairs
    return snp_df


def per_snp_summary_for_region(mat, annotated_pairs, snp_labels, region) -> pd.DataFrame:
    mask = (annotated_pairs["region1"] == region) & (annotated_pairs["region2"] == region)
    row_idx = annotated_pairs.index[mask].values
    n_within = len(row_idx)
    if n_within == 0:
        return pd.DataFrame()
    counts = np.asarray(mat[row_idx, :].sum(axis=0)).ravel()
    snp_df = parse_snp_labels(snp_labels)
    snp_df["region"] = region
    snp_df["n_pairs_ibd"] = counts
    snp_df["n_pairs_total"] = n_within
    snp_df["frac_pairs_ibd"] = counts / n_within
    return snp_df


def per_snp_summary_between_regions(mat, annotated_pairs, snp_labels, region_a, region_b) -> pd.DataFrame:
    if region_a == region_b:
        mask = (annotated_pairs["region1"] == region_a) & (annotated_pairs["region2"] == region_a)
    else:
        mask = (
            ((annotated_pairs["region1"] == region_a) & (annotated_pairs["region2"] == region_b)) |
            ((annotated_pairs["region1"] == region_b) & (annotated_pairs["region2"] == region_a))
        )
    row_idx = annotated_pairs.index[mask].values
    n_pairs = len(row_idx)
    if n_pairs == 0:
        return pd.DataFrame()
    counts = np.asarray(mat[row_idx, :].sum(axis=0)).ravel()
    snp_df = parse_snp_labels(snp_labels)
    snp_df["region_a"] = region_a
    snp_df["region_b"] = region_b
    snp_df["n_pairs_ibd"] = counts
    snp_df["n_pairs_total"] = n_pairs
    snp_df["frac_pairs_ibd"] = counts / n_pairs
    return snp_df


def annotate_pairs_with_regions(pair_summary, meta, region_col="region") -> pd.DataFrame:
    meta_idx = meta.set_index("sample")[region_col].to_dict()
    pair_summary = pair_summary.copy()
    pair_summary["region1"] = pair_summary["sample1"].map(meta_idx).fillna("unknown")
    pair_summary["region2"] = pair_summary["sample2"].map(meta_idx).fillna("unknown")
    pair_summary["same_region"] = pair_summary["region1"] == pair_summary["region2"]
    return pair_summary


def within_between_region_ibd(annotated_pairs) -> pd.DataFrame:
    return (
        annotated_pairs.groupby("same_region")["frac_ibd"]
        .agg(["mean", "median", "std", "count"])
        .rename(index={True: "within_region", False: "between_region"})
        .reset_index()
        .rename(columns={"same_region": "comparison_type"})
    )


def pairwise_region_ibd(annotated_pairs) -> pd.DataFrame:
    df = annotated_pairs.copy()
    df["reg_a"] = np.where(df["region1"] <= df["region2"], df["region1"], df["region2"])
    df["reg_b"] = np.where(df["region1"] <= df["region2"], df["region2"], df["region1"])
    return (
        df.groupby(["reg_a", "reg_b"])["frac_ibd"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
    )


def per_chr_ibd(mat, snp_labels: list) -> pd.DataFrame:
    snp_df = parse_snp_labels(snp_labels)
    records = []
    for chrom, grp in snp_df.groupby("chr"):
        col_idx = grp.index.values
        sub = mat[:, col_idx]
        n_snps_chr = len(col_idx)
        frac_per_pair = np.asarray(sub.sum(axis=1)).ravel() / n_snps_chr
        records.append({
            "chr": chrom, "n_snps": n_snps_chr,
            "mean_frac_ibd": frac_per_pair.mean(),
            "median_frac_ibd": np.median(frac_per_pair),
            "std_frac_ibd": frac_per_pair.std(),
        })
    return pd.DataFrame(records).sort_values("chr")
