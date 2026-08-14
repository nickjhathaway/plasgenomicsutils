#!/usr/bin/env python
"""Per-pair IBD fraction and SNP density from hmmibd-rs blocks."""

from __future__ import annotations

import argparse
from pathlib import Path

from ...lib import ibd_fraction as F
from ...lib.reference import DEFAULT_REFERENCE, available_references, get_reference
from ...lib.ibd_matrix import IBD_MIN_BLOCK_KB, IBD_MIN_BLOCK_SNP
from ...lib.vcf_io import positions_frame
from ...utils.small_utils import Utils


def get_parser_fraction_and_snp_density() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils ibd_fraction_and_snp_density",
        description="Per-pair IBD fraction (callable denominator) and SNP density",
    )
    p.add_argument("--blocks", required=True, help="hmmibd-rs blocks TSV")
    p.add_argument("--snps", required=True,
                   help="VCF/BED of the SNP set (needed for callable denominator)")
    p.add_argument("--snp-format", choices=["vcf", "bed"], default="vcf")
    p.add_argument("--reference", default=DEFAULT_REFERENCE,
                   help=f"Reference genome for chromosome lengths (default: {DEFAULT_REFERENCE}; "
                        f"available: {', '.join(available_references())})")
    p.add_argument("--bp-per-cm", type=float, default=None,
                   help="Constant genetic-map rate (bp/cM). Default: the reference's rate "
                        "(Pf3D7 = 15000).")
    p.add_argument("--min-snp", type=int, default=15, help="SNP-per-1cM-window floor to report")
    p.add_argument("--sep", default="\t")
    p.add_argument("--output", default="ibd_frac", help="Output prefix")
    p.add_argument("--min-block-snp", type=int, default=IBD_MIN_BLOCK_SNP,
                   help="Drop IBD segments with fewer than this many SNPs; short, SNP-poor "
                        "segments are commonly spurious (default: %(default)s, 0 disables)")
    p.add_argument("--min-block-kb", type=float, default=IBD_MIN_BLOCK_KB,
                   help="Drop IBD segments shorter than this many kb "
                        "(default: %(default)s, 0 disables)")
    return p


def parse_args_fraction_and_snp_density():
    return get_parser_fraction_and_snp_density().parse_args()


def fraction_and_snp_density():
    args = parse_args_fraction_and_snp_density()
    ref = get_reference(args.reference)
    bp_per_cm = args.bp_per_cm if args.bp_per_cm is not None else ref.bp_per_cm

    pos_df = positions_frame(args.snps, args.snp_format)
    span_df = F.callable_spans(pos_df, ref)
    callable_bp = int(span_df["span_bp"].sum())
    full_bp = int(span_df["full_bp"].sum())
    callable_cm = callable_bp / bp_per_cm
    full_cm = full_bp / bp_per_cm

    print("--- Callable genome (first->last SNP per chromosome) ---")
    print(span_df.to_string(index=False))
    print(f"\nCallable: {callable_bp:,} bp = {callable_cm:,.1f} cM"
          f"   ({100*callable_bp/full_bp:.1f}% of full {full_bp:,} bp)")

    print("\n--- Per-pair IBD fraction (callable denominator) ---")
    pair_df = F.per_pair_fraction(args.blocks, args.sep, bp_per_cm, callable_cm,
                                  min_block_snp=args.min_block_snp,
                                  min_block_kb=args.min_block_kb)
    pair_df["ibd_fraction_full_genome"] = pair_df["total_ibd_cm"] / full_cm
    # sample1/sample2 first: the table is an edge list before it is anything else
    lead = ["pair", "sample1", "sample2"]
    pair_df = pair_df[lead + [c for c in pair_df.columns if c not in lead]]
    n_with = int((pair_df["total_ibd_bp"] > 0).sum())
    print(f"  {len(pair_df):,} total pairs;  {n_with:,} with IBD;  "
          f"{len(pair_df) - n_with:,} zero-IBD")
    print(pair_df["ibd_fraction_accessible"].describe().to_string())
    out_pairs = f"{args.output}.pair_ibd_fraction.tsv.gz"
    Utils.write_tsv_gz(pair_df, out_pairs)
    print(f"  -> {out_pairs}")

    print("\n--- SNP density ---")
    summary, win_df = F.snp_density(pos_df, ref, bp_per_cm, callable_cm, full_cm, args.min_snp)
    for k, v in summary.items():
        print(f"  {k}: {v}")
    Path(f"{args.output}.snp_density.txt").write_text(
        "\n".join(f"{k}\t{v}" for k, v in summary.items()) + "\n")
    out_win = f"{args.output}.snp_density_windows.tsv.gz"
    Utils.write_tsv_gz(win_df, out_win)
    print(f"  -> {args.output}.snp_density.txt")
    print(f"  -> {out_win}")


if __name__ == "__main__":
    fraction_and_snp_density()
