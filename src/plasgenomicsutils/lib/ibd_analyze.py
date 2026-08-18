"""Downstream analysis of the binary IBD matrix.

Per-pair / per-SNP / per-group / per-chromosome IBD summaries.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .ibd_selection import parse_snp_labels


def parse_pair_labels(pair_labels: list) -> pd.DataFrame:
    """Split ``sample1__sample2`` pair labels into a frame of ``pair, sample1, sample2``."""
    rows = []
    for label in pair_labels:
        s1, s2 = label.split("__", 1)
        rows.append({"pair": label, "sample1": s1, "sample2": s2})
    return pd.DataFrame(rows)


def per_pair_summary(mat, pair_labels: list) -> pd.DataFrame:
    """One row per sample pair: how many SNPs it is IBD at, and what fraction of the panel.

    ``frac_ibd`` divides by the number of SNPs in the matrix, so it is a property of the
    panel as much as of the pair -- comparable within a run, not between panels. For a
    length-based fraction use ``ibd_fraction_and_snp_density`` instead, which divides by
    callable cM.
    """
    n_snps = mat.shape[1]
    counts = np.asarray(mat.sum(axis=1)).ravel()
    pair_df = parse_pair_labels(pair_labels)
    pair_df["n_ibd_snps"] = counts
    pair_df["frac_ibd"] = counts / n_snps
    return pair_df


def per_snp_summary(mat, snp_labels: list) -> pd.DataFrame:
    """One row per SNP: how many pairs are IBD there, over every pair in the matrix."""
    n_pairs = mat.shape[0]
    counts = np.asarray(mat.sum(axis=0)).ravel()
    snp_df = parse_snp_labels(snp_labels)
    snp_df["n_pairs_ibd"] = counts
    snp_df["frac_pairs_ibd"] = counts / n_pairs
    return snp_df


def per_snp_summary_for_group(mat, annotated_pairs, snp_labels, group) -> pd.DataFrame:
    """Per-SNP IBD within one group: only pairs whose members are both in ``group``.

    Returns an empty frame when the group has no within-group pair, which is the honest
    answer for a group of one rather than a row of zeros.
    """
    mask = (annotated_pairs["group1"] == group) & (annotated_pairs["group2"] == group)
    row_idx = annotated_pairs.index[mask].values
    n_within = len(row_idx)
    if n_within == 0:
        return pd.DataFrame()
    counts = np.asarray(mat[row_idx, :].sum(axis=0)).ravel()
    snp_df = parse_snp_labels(snp_labels)
    snp_df["group"] = group
    snp_df["n_pairs_ibd"] = counts
    snp_df["n_pairs_total"] = n_within
    snp_df["frac_pairs_ibd"] = counts / n_within
    return snp_df


def per_snp_summary_between_groups(mat, annotated_pairs, snp_labels, group_a, group_b) -> pd.DataFrame:
    """Per-SNP IBD for one pair of groups, in either order.

    ``group_a == group_b`` gives that group's within-group pairs, so the same function
    fills the diagonal and the off-diagonal of a group x group panel.
    """
    if group_a == group_b:
        mask = (annotated_pairs["group1"] == group_a) & (annotated_pairs["group2"] == group_a)
    else:
        mask = (
            ((annotated_pairs["group1"] == group_a) & (annotated_pairs["group2"] == group_b)) |
            ((annotated_pairs["group1"] == group_b) & (annotated_pairs["group2"] == group_a))
        )
    row_idx = annotated_pairs.index[mask].values
    n_pairs = len(row_idx)
    if n_pairs == 0:
        return pd.DataFrame()
    counts = np.asarray(mat[row_idx, :].sum(axis=0)).ravel()
    snp_df = parse_snp_labels(snp_labels)
    snp_df["group_a"] = group_a
    snp_df["group_b"] = group_b
    snp_df["n_pairs_ibd"] = counts
    snp_df["n_pairs_total"] = n_pairs
    snp_df["frac_pairs_ibd"] = counts / n_pairs
    return snp_df


def annotate_pairs_with_groups(pair_summary, meta, group_col="group") -> pd.DataFrame:
    """Attach each pair's two group labels, and whether they match.

    Samples absent from ``meta`` are labelled ``"unknown"`` rather than dropped: a pair
    silently disappearing would shrink the denominator of every fraction computed from
    this frame.
    """
    meta_idx = meta.set_index("sample")[group_col].to_dict()
    pair_summary = pair_summary.copy()
    pair_summary["group1"] = pair_summary["sample1"].map(meta_idx).fillna("unknown")
    pair_summary["group2"] = pair_summary["sample2"].map(meta_idx).fillna("unknown")
    pair_summary["same_group"] = pair_summary["group1"] == pair_summary["group2"]
    return pair_summary


def within_between_group_ibd(annotated_pairs) -> pd.DataFrame:
    """Mean / median / sd of pair sharing, split into within-group and between-group."""
    return (
        annotated_pairs.groupby("same_group")["frac_ibd"]
        .agg(["mean", "median", "std", "count"])
        .rename(index={True: "within_group", False: "between_group"})
        .reset_index()
        .rename(columns={"same_group": "comparison_type"})
    )


def pairwise_group_ibd(annotated_pairs) -> pd.DataFrame:
    """Pair sharing summarised for every group pair, each unordered pair once.

    The two group labels are sorted per row before grouping, so A-vs-B and B-vs-A land on
    the same row instead of being counted as two.
    """
    df = annotated_pairs.copy()
    df["reg_a"] = np.where(df["group1"] <= df["group2"], df["group1"], df["group2"])
    df["reg_b"] = np.where(df["group1"] <= df["group2"], df["group2"], df["group1"])
    return (
        df.groupby(["reg_a", "reg_b"])["frac_ibd"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
    )


def per_chr_ibd(mat, snp_labels: list) -> pd.DataFrame:
    """Pair sharing summarised per chromosome, with the SNP count behind each row.

    ``n_snps`` matters here: chromosomes carry very different numbers of SNPs, so a
    chromosome's mean is estimated far more precisely on some than on others.
    """
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
