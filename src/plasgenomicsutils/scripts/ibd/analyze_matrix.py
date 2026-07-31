#!/usr/bin/env python
"""Downstream per-pair / per-SNP / per-region / per-chromosome IBD summaries."""

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
    p.add_argument("--matrix", required=True, help="Prefix from build_ibd_matrix")
    p.add_argument("--meta", help="Sample metadata CSV (columns: sample, <region_col>)")
    p.add_argument("--region-col", default="region",
                   help="Column in metadata to use as region (default: region)")
    p.add_argument("--output", default="ibd_analysis", help="Output prefix")
    p.add_argument("--pairwise-region-snp", action="store_true",
                   help="Also compute per-SNP IBD for every pairwise region combination. "
                        "Output: <output>.per_snp_pairwise_region.tsv.gz. Requires --meta.")
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
        print("\n(No --meta provided; skipping region-level analysis)")
        return

    print(f"\n--- Region-level IBD (region_col='{args.region_col}') ---")
    meta = pd.read_csv(args.meta)
    annotated = A.annotate_pairs_with_regions(pair_sum, meta, args.region_col)
    _save(annotated, f"{args.output}.per_pair_annotated.tsv.gz")

    wb = A.within_between_region_ibd(annotated)
    print(wb.to_string(index=False))
    _save(wb, f"{args.output}.within_between_region.tsv.gz")

    _save(A.pairwise_region_ibd(annotated), f"{args.output}.pairwise_region_ibd.tsv.gz")

    print("\n--- Per-region per-SNP IBD frequency ---")
    regions = sorted(set(annotated["region1"].tolist() + annotated["region2"].tolist()) - {"unknown"})
    print(f"Found {len(regions)} regions: {', '.join(regions)}")
    region_snp_dfs = []
    for region in regions:
        region_snp = A.per_snp_summary_for_region(mat, annotated, snp_labels, region)
        if region_snp.empty:
            print(f"  Skipping '{region}' — no within-region pairs")
            continue
        print(f"  Processed '{region}' ({region_snp['n_pairs_total'].iloc[0]:,} pairs)")
        region_snp_dfs.append(region_snp)
    if region_snp_dfs:
        _save(pd.concat(region_snp_dfs, ignore_index=True), f"{args.output}.per_snp_per_region.tsv.gz")

    if args.pairwise_region_snp:
        print("\n--- Pairwise region per-SNP IBD frequency ---")
        region_pairs = [(ra, rb) for i, ra in enumerate(regions) for rb in regions[i:]]
        print(f"  {len(region_pairs)} region pairs")
        pairwise_snp_dfs = []
        for ra, rb in region_pairs:
            df = A.per_snp_summary_between_regions(mat, annotated, snp_labels, ra, rb)
            if df.empty:
                print(f"  Skipping '{ra}' x '{rb}' — no pairs")
                continue
            print(f"  '{ra}' x '{rb}': {df['n_pairs_total'].iloc[0]:,} pairs")
            pairwise_snp_dfs.append(df)
        if pairwise_snp_dfs:
            _save(pd.concat(pairwise_snp_dfs, ignore_index=True),
                  f"{args.output}.per_snp_pairwise_region.tsv.gz")

    print(f"\nAll outputs written with prefix '{args.output}.*'")


if __name__ == "__main__":
    analyze_matrix()
