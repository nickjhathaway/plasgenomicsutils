#!/usr/bin/env python
"""IBD-based selection test statistic (XiR,s), genome-wide and per-region."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ...lib import ibd_matrix
from ...lib import ibd_selection as S
from ...utils.small_utils import Utils


def get_parser_selection_statistic() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils ibd_selection_statistic",
        description="Compute IBD-based selection test statistic (XiR,s)",
    )
    p.add_argument("--matrix", required=True, help="Prefix from build_ibd_matrix")
    p.add_argument("--af", required=True,
                   help="TSV (snp_id, af): global allele frequencies (required; "
                        "generate with compute_allele_freqs --zero-based). AFs are not "
                        "estimated from the matrix.")
    p.add_argument("--af-region", default=None,
                   help="TSV (region, snp_id, af): per-region AFs; falls back to --af.")
    p.add_argument("--meta", default=None,
                   help="Sample metadata CSV with columns: sample, <region-col>")
    p.add_argument("--region-col", default="region",
                   help="Column in --meta to use as region (default: region)")
    p.add_argument("--n-bins", type=int, default=100,
                   help="Number of AF bins for normalisation (default: 100)")
    p.add_argument("--alpha", type=float, default=0.05,
                   help="Genome-wide significance level (default: 0.05)")
    p.add_argument("--output", default="ibd_selection", help="Output prefix")
    return p


def parse_args_selection_statistic():
    return get_parser_selection_statistic().parse_args()


def _save(df, path):
    Utils.write_tsv_gz(df, path)
    print(f"  -> {path}")


def selection_statistic():
    args = parse_args_selection_statistic()

    mat, pair_labels, snp_labels = ibd_matrix.load_matrix(args.matrix)
    print(f"Loaded matrix: {mat.shape[0]:,} pairs x {mat.shape[1]:,} SNPs")
    snp_df = S.parse_snp_labels(snp_labels)

    print("\n--- Allele frequencies ---")
    global_af = S.load_global_af(args.af, snp_labels)
    region_af_table = S.load_region_af_table(args.af_region) if args.af_region else None

    print("\n--- Global selection statistic ---")
    stats, bin_df = S.compute_selection_statistic(mat, global_af, n_bins=args.n_bins, label="global")
    out_df, threshold = S.assemble_output(snp_df, stats, args.alpha)
    n_valid = int((~out_df["neg_log10_p"].isna()).sum())
    n_sig = int(out_df["significant"].sum())
    print(f"  Valid SNPs: {n_valid:,}  |  Bonferroni threshold: {threshold:.3f}  |  Significant: {n_sig:,}")

    _save(out_df, f"{args.output}.global.selection_stats.tsv.gz")
    _save(bin_df, f"{args.output}.global.bin_stats.tsv.gz")
    Path(f"{args.output}.global.threshold.txt").write_text(
        "region\talpha\tn_tests\tneg_log10_p_threshold\n"
        f"global\t{args.alpha}\t{n_valid}\t{threshold:.6f}\n"
    )
    print(f"  -> {args.output}.global.threshold.txt")

    if not args.meta:
        print("\n(No --meta provided; skipping per-region analysis)")
        return

    print(f"\n--- Per-region selection statistic (region_col='{args.region_col}') ---")
    meta = pd.read_csv(args.meta)
    regions = sorted(meta[args.region_col].dropna().unique())
    print(f"Found {len(regions)} regions: {', '.join(regions)}")

    region_stats_dfs, region_bin_dfs, threshold_rows = [], [], []
    for region in regions:
        print(f"\n  Region: {region}")
        row_idx = S.within_region_row_indices(pair_labels, meta, args.region_col, region)
        if len(row_idx) < 2:
            print(f"    Skipping — only {len(row_idx)} within-region pair(s)")
            continue
        mat_sub = mat[row_idx, :]
        af_reg = S.get_af_for_region(region, snp_labels, region_af_table, global_af)
        stats_r, bin_df_r = S.compute_selection_statistic(mat_sub, af_reg, n_bins=args.n_bins, label=region)
        out_r, thresh_r = S.assemble_output(snp_df, stats_r, args.alpha, region=region)
        n_valid_r = int((~out_r["neg_log10_p"].isna()).sum())
        n_sig_r = int(out_r["significant"].sum())
        print(f"    Valid SNPs: {n_valid_r:,}  |  Bonferroni threshold: {thresh_r:.3f}  |  Significant: {n_sig_r:,}")
        bin_df_r.insert(0, "region", region)
        region_stats_dfs.append(out_r)
        region_bin_dfs.append(bin_df_r)
        threshold_rows.append({
            "region": region, "alpha": args.alpha,
            "n_tests": n_valid_r, "neg_log10_p_threshold": thresh_r,
        })

    if region_stats_dfs:
        print("\n--- Saving per-region outputs ---")
        _save(pd.concat(region_stats_dfs, ignore_index=True), f"{args.output}.per_region.selection_stats.tsv.gz")
        _save(pd.concat(region_bin_dfs, ignore_index=True), f"{args.output}.per_region.bin_stats.tsv.gz")
        thresh_path = f"{args.output}.per_region.threshold.txt"
        pd.DataFrame(threshold_rows).to_csv(thresh_path, sep="\t", index=False)
        print(f"  -> {thresh_path}")
    else:
        print("  No regions had sufficient pairs; no per-region output written.")


if __name__ == "__main__":
    selection_statistic()
