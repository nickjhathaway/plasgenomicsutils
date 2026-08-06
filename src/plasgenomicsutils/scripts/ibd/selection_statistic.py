#!/usr/bin/env python
"""IBD-based selection test statistic (XiR,s), genome-wide and per-group."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ...lib import ibd_matrix
from ...lib import ibd_selection as S
from ...utils.small_utils import Utils


def get_parser_selection_statistic() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils ibd_selection_statistic",
        description="Compute IBD-based selection test statistic (XiR,s)",
    )
    p.add_argument("--matrix", required=True,
                   help="Matrix from build_ibd_matrix: the prefix or the .npz path")
    p.add_argument("--af", required=True,
                   help="TSV (snp_id, af): global allele frequencies (required; generate "
                        "with compute_allele_freqs). AFs are not estimated from the matrix.")
    p.add_argument("--af-group", default=None,
                   help="TSV (group, snp_id, af): per-group AFs; falls back to --af.")
    p.add_argument("--meta", default=None,
                   help="Sample metadata CSV with columns: sample, <group-col>")
    p.add_argument("--group-col", default="group",
                   help="Column in --meta to use as group (default: group)")
    p.add_argument("--xirs-variant", choices=("corrected", "published"), default="corrected",
                   help="'published' reproduces the recipe in Henden et al. and its "
                        "implementations (isoRelate `iRfunction`, ibdutils `calc_xirs`): it "
                        "centres each SNP and then sums that same SNP, which cancels to zero, "
                        "so the output is floating-point residue and two implementations of "
                        "it disagree. Use only to reproduce prior work. 'corrected' (default) "
                        "omits the per-SNP centring, so the statistic measures what it "
                        "claims: excess IBD sharing at a locus (default: %(default)s)")
    p.add_argument("--tail", choices=("upper", "two-sided"), default="upper",
                   help="'upper' asks only whether a locus is shared MORE than expected, "
                        "which is what a positive-selection scan means. 'two-sided' scores a "
                        "sharing deficit identically to an excess (what z^2 -> chi2(1) does, "
                        "and what the published recipe uses) (default: %(default)s)")
    p.add_argument("--n-bins", type=int, default=100,
                   help="Number of AF bins for normalisation (default: 100)")
    p.add_argument("--alpha", type=float, default=0.05,
                   help="Family-wise level for the Bonferroni threshold, which controls "
                        "the chance of even one false positive (default: %(default)s)")
    p.add_argument("--fdr-alpha", type=float, default=None,
                   help="q-value cutoff for the Benjamini-Hochberg column, which instead "
                        "controls the expected share of false positives among the SNPs "
                        "called. Both are always written; this only sets where the "
                        "`significant_fdr` flag falls (default: same as --alpha)")
    p.add_argument("--permute", type=int, default=0, metavar="N",
                   help="Build a null by sliding every pair's IBD segments to random "
                        "positions N times and recomputing. Unlike Bonferroni and BH this "
                        "assumes nothing about the chi2 fit and accounts for one segment "
                        "spanning many SNPs. Yields a family-wise threshold (the 1-alpha "
                        "quantile of the genome-wide maxima), per-SNP empirical p-values "
                        "and a Benjamini-Hochberg pass over those -- the only FDR here "
                        "resting on calibrated p-values. Expect all of it to be far "
                        "stricter. 200 is enough for alpha=0.05 (default: %(default)s, off)")
    p.add_argument("--empirical-pool", choices=("global", "bin"), default="global",
                   help="Which null --permute pools to get each SNP's empirical p-value. "
                        "'global' uses every SNP's nulls, which resolves ~n_bins times "
                        "finer and is what a genome-wide FDR needs, but assumes the null "
                        "has the same shape in every MAF bin. 'bin' stays inside the "
                        "SNP's own bin, dropping that assumption for a much coarser "
                        "p-value. Both are always written (`p_empirical`, "
                        "`p_empirical_binned`); this picks the one `q_empirical` uses "
                        "(default: %(default)s)")
    p.add_argument("--permute-seed", type=int, default=0,
                   help="RNG seed for --permute (default: %(default)s)")
    p.add_argument("--output", default="ibd_selection", help="Output prefix")
    return p


def parse_args_selection_statistic():
    return get_parser_selection_statistic().parse_args()


def _permute(args, mat, af, indent="  "):
    """Run the permutation for one scan, or return None when --permute is off."""
    if not args.permute:
        return None
    print(f"{indent}permuting ({args.permute} replicates)...", flush=True)
    return S.permutation_null(mat, af, n_perm=args.permute, n_bins=args.n_bins,
                              alpha=args.alpha, seed=args.permute_seed,
                              variant=args.xirs_variant, tail=args.tail)


def _report(info, label, indent="  "):
    """One line per scan: the corrections, and the calibration behind them."""
    lam = info["lambda_gc"]
    if info.get("xirs_variant") == "published":
        print(f"{indent}WARNING: --xirs-variant published sums each SNP after centring it, "
              f"which cancels to zero, so these numbers are floating-point residue rather "
              f"than a statistic. For reproducing prior work only.")
    print(f"{indent}Valid SNPs: {info['n_tests']:,}"
          f"  |  Bonferroni (a={info['alpha']}) >= {info['neg_log10_p_threshold']:.2f}: "
          f"{info['n_significant']:,}"
          f"  |  BH FDR (q<{info['fdr_alpha']}): {info['n_significant_fdr']:,}")
    if np.isfinite(info.get("neg_log10_p_perm_threshold", np.nan)):
        print(f"{indent}Permutation ({info['n_perm']} replicates), built from the data so "
              f"no chi2(1) assumption -- prefer these two where they disagree with the above:")
        print(f"{indent}  FWER {info['alpha']} >= {info['neg_log10_p_perm_threshold']:.2f}: "
              f"{info['n_significant_perm']:,}")
        print(f"{indent}  BH on empirical p (q<{info['fdr_alpha']}): "
              f"{info['n_significant_fdr_perm']:,}"
              + (f" (>= {info['neg_log10_p_emp_fdr_threshold']:.2f})"
                 if np.isfinite(info["neg_log10_p_emp_fdr_threshold"]) else ""))
        qf, dead = info["q_empirical_floor"], info["frac_q_unreachable"]
        if np.isfinite(qf) and qf > info["fdr_alpha"] / 10:
            need = int(np.ceil(info["n_perm"] * qf / (info["fdr_alpha"] / 10) / 100) * 100)
            k = int(np.ceil(qf / info["fdr_alpha"]))
            print(f"{indent}  warning: with {info['n_perm']} replicates a lone top SNP "
                  f"bottoms out at q={qf:.3f}, so nothing clears q<{info['fdr_alpha']} "
                  f"unless at least {k} SNPs tie at the resolution limit (BH divides by "
                  f"rank). --permute {need} would give an order of magnitude of headroom.")
        if np.isfinite(dead) and dead > 0.01:
            print(f"{indent}  warning: {dead:.0%} of SNPs sit in a MAF bin too small to "
                  f"reach q<{info['fdr_alpha']} however extreme their score, so the "
                  f"empirical FDR is blind to them. Raise --permute, or lower --n-bins "
                  f"to put more SNPs behind each p-value.")
        lo, hi = info["perm_bin_tail_min"], info["perm_bin_tail_max"]
        if np.isfinite(hi) and not (0.005 <= lo and hi <= 0.02):
            off = (max(hi / 0.01, 0.01 / lo) if lo > 0 else np.inf)
            size = f"~{off:.1f}x" if np.isfinite(off) else "an unbounded factor"
            print(f"{indent}  note: MAF bins' null tails run {lo:.3f}-{hi:.3f} where "
                  f"pooling wants 0.010, so `p_empirical` is off by up to {size} for "
                  f"SNPs in the outlying bins -- too small in the heavy-tailed ones. Top "
                  f"hits are unaffected at that size; SNPs near the q cutoff are not. "
                  f"`p_empirical_binned` needs no pooling; "
                  + ("it already drives `q_empirical`"
                     if info["empirical_pool"] == "bin"
                     else "compare against it, or switch with --empirical-pool bin"))
    if 0 < info["n_bins_used"] and info["largest_bin_frac"] > 0.1:
        print(f"{indent}MAF ties collapsed the binning to {info['n_bins_used']} non-empty "
              f"bins, the largest holding {info['largest_bin_frac']:.0%} of SNPs -- fewer "
              f"than requested, which coarsens the within-bin standardisation the "
              f"statistic itself rests on, not just the p-values.")
    note = ""
    if np.isfinite(lam) and not (0.8 <= lam <= 1.25):
        note = ("  <- far from 1: the chi2(1) null does not fit, so Bonferroni and BH both "
                "rest on miscalibrated p-values")
        note += ("; use the permutation threshold" if np.isfinite(
            info.get("neg_log10_p_perm_threshold", np.nan))
            else "; re-run with --permute, or rank SNPs and merge peaks instead")
    print(f"{indent}Genomic inflation lambda = {lam:.2f}{note}")


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
    group_af_table = S.load_group_af_table(args.af_group) if args.af_group else None

    print("\n--- Global selection statistic ---")
    stats, bin_df = S.compute_selection_statistic(mat, global_af, n_bins=args.n_bins,
                                                 label="global", variant=args.xirs_variant,
                                                 tail=args.tail)
    perm = _permute(args, mat, global_af)
    out_df, info = S.assemble_output(snp_df, stats, args.alpha, fdr_alpha=args.fdr_alpha,
                                     perm=perm, pool=args.empirical_pool,
                                     variant=args.xirs_variant, tail=args.tail)
    n_valid = info["n_tests"]
    _report(info, "global")

    _save(out_df, f"{args.output}.global.selection_stats.tsv.gz")
    _save(bin_df, f"{args.output}.global.bin_stats.tsv.gz")
    pd.DataFrame([dict(group="global", **info)]).to_csv(
        f"{args.output}.global.threshold.txt", sep="\t", index=False)
    print(f"  -> {args.output}.global.threshold.txt")

    if not args.meta:
        print("\n(No --meta provided; skipping per-group analysis)")
        return

    print(f"\n--- Per-group selection statistic (group_col='{args.group_col}') ---")
    meta = Utils.read_table(args.meta)   # auto-detect tab / comma
    groups = sorted(meta[args.group_col].dropna().unique())
    print(f"Found {len(groups)} groups: {', '.join(groups)}")

    group_stats_dfs, group_bin_dfs, threshold_rows = [], [], []
    for group in groups:
        print(f"\n  Group: {group}")
        row_idx = S.within_group_row_indices(pair_labels, meta, args.group_col, group)
        if len(row_idx) < 2:
            print(f"    Skipping — only {len(row_idx)} within-group pair(s)")
            continue
        mat_sub = mat[row_idx, :]
        af_reg = S.get_af_for_group(group, snp_labels, group_af_table, global_af)
        stats_r, bin_df_r = S.compute_selection_statistic(mat_sub, af_reg, n_bins=args.n_bins,
                                                         label=group, variant=args.xirs_variant,
                                                         tail=args.tail)
        perm_r = _permute(args, mat_sub, af_reg, indent="    ")
        out_r, info_r = S.assemble_output(snp_df, stats_r, args.alpha, group=group,
                                          fdr_alpha=args.fdr_alpha, perm=perm_r,
                                          pool=args.empirical_pool,
                                          variant=args.xirs_variant, tail=args.tail)
        _report(info_r, group, indent="    ")
        bin_df_r.insert(0, "group", group)
        group_stats_dfs.append(out_r)
        group_bin_dfs.append(bin_df_r)
        threshold_rows.append(dict(group=group, **info_r))

    if group_stats_dfs:
        print("\n--- Saving per-group outputs ---")
        _save(pd.concat(group_stats_dfs, ignore_index=True), f"{args.output}.per_group.selection_stats.tsv.gz")
        _save(pd.concat(group_bin_dfs, ignore_index=True), f"{args.output}.per_group.bin_stats.tsv.gz")
        thresh_path = f"{args.output}.per_group.threshold.txt"
        pd.DataFrame(threshold_rows).to_csv(thresh_path, sep="\t", index=False)
        print(f"  -> {thresh_path}")
    else:
        print("  No groups had sufficient pairs; no per-group output written.")


if __name__ == "__main__":
    selection_statistic()
