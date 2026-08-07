#!/usr/bin/env python
"""Downstream per-pair / per-SNP / per-group / per-chromosome IBD summaries."""

from __future__ import annotations

import argparse

import pandas as pd

from ...lib import ibd_analyze as A
from ...lib import ibd_matrix
from ...utils.small_utils import Utils


def get_parser_analyze_matrix() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils analyze_ibd_matrix",
        description="Downstream analysis of the binary IBD matrix",
    )
    p.add_argument("--matrix", required=True,
                   help="Matrix from build_ibd_matrix: the prefix or the .npz path")
    p.add_argument("--meta", help="Sample metadata CSV (columns: sample, <group_col>)")
    p.add_argument("--group-col", default="group",
                   help="Column in metadata to use as group (default: group)")
    p.add_argument("--output", default="ibd_analysis", help="Output prefix")
    p.add_argument("--pairwise-group-snp", action="store_true",
                   help="Also compute per-SNP IBD for every pairwise group combination. "
                        "Output: <output>.per_snp_pairwise_group.tsv.gz. Requires --meta.")
    return p


def parse_args_analyze_matrix():
    return get_parser_analyze_matrix().parse_args()


def _save(df, path):
    Utils.write_tsv_gz(df, path)
    print(f"  -> {path}")


def analyze_matrix():
    args = parse_args_analyze_matrix()

    mat, pair_labels, snp_labels = ibd_matrix.load_matrix(args.matrix)
    print(f"Loaded matrix: {mat.shape[0]:,} pairs x {mat.shape[1]:,} SNPs")

    print("\n--- Per-pair IBD summary ---")
    pair_sum = A.per_pair_summary(mat, pair_labels)
    print(pair_sum["frac_ibd"].describe().to_string())
    _save(pair_sum, f"{args.output}.per_pair.tsv.gz")

    print("\n--- Per-SNP IBD frequency (global) ---")
    snp_sum = A.per_snp_summary(mat, snp_labels)
    print(snp_sum["frac_pairs_ibd"].describe().to_string())
    _save(snp_sum, f"{args.output}.per_snp.tsv.gz")

    print("\n--- Per-chromosome mean IBD ---")
    chr_sum = A.per_chr_ibd(mat, snp_labels)
    print(chr_sum.to_string(index=False))
    _save(chr_sum, f"{args.output}.per_chr.tsv.gz")

    if not args.meta:
        print("\n(No --meta provided; skipping group-level analysis)")
        return

    print(f"\n--- Group-level IBD (group_col='{args.group_col}') ---")
    meta = Utils.read_table(args.meta)   # auto-detect tab / comma
    annotated = A.annotate_pairs_with_groups(pair_sum, meta, args.group_col)
    _save(annotated, f"{args.output}.per_pair_annotated.tsv.gz")

    wb = A.within_between_group_ibd(annotated)
    print(wb.to_string(index=False))
    _save(wb, f"{args.output}.within_between_group.tsv.gz")

    _save(A.pairwise_group_ibd(annotated), f"{args.output}.pairwise_group_ibd.tsv.gz")

    print("\n--- Per-group per-SNP IBD frequency ---")
    groups = sorted(set(annotated["group1"].tolist() + annotated["group2"].tolist()) - {"unknown"})
    print(f"Found {len(groups)} groups: {', '.join(groups)}")
    group_snp_dfs = []
    for group in groups:
        group_snp = A.per_snp_summary_for_group(mat, annotated, snp_labels, group)
        if group_snp.empty:
            print(f"  Skipping '{group}' — no within-group pairs")
            continue
        print(f"  Processed '{group}' ({group_snp['n_pairs_total'].iloc[0]:,} pairs)")
        group_snp_dfs.append(group_snp)
    if group_snp_dfs:
        _save(pd.concat(group_snp_dfs, ignore_index=True), f"{args.output}.per_snp_per_group.tsv.gz")

    if args.pairwise_group_snp:
        print("\n--- Pairwise group per-SNP IBD frequency ---")
        group_pairs = [(ra, rb) for i, ra in enumerate(groups) for rb in groups[i:]]
        print(f"  {len(group_pairs)} group pairs")
        pairwise_snp_dfs = []
        for ra, rb in group_pairs:
            df = A.per_snp_summary_between_groups(mat, annotated, snp_labels, ra, rb)
            if df.empty:
                print(f"  Skipping '{ra}' x '{rb}' — no pairs")
                continue
            print(f"  '{ra}' x '{rb}': {df['n_pairs_total'].iloc[0]:,} pairs")
            pairwise_snp_dfs.append(df)
        if pairwise_snp_dfs:
            _save(pd.concat(pairwise_snp_dfs, ignore_index=True),
                  f"{args.output}.per_snp_pairwise_group.tsv.gz")

    print(f"\nAll outputs written with prefix '{args.output}.*'")


if __name__ == "__main__":
    analyze_matrix()
