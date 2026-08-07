#!/usr/bin/env python
"""Build a binary (pairs x SNPs) IBD matrix from hmmibd-rs blocks + a SNP panel."""

from __future__ import annotations

import argparse

from ...lib import ibd_matrix
from ...lib.ibd_matrix import IBD_MIN_BLOCK_KB, IBD_MIN_BLOCK_SNP
from ...lib.vcf_io import SnpPanel


def get_parser_build_matrix() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils build_ibd_matrix",
        description="Build binary IBD matrix from hmmibd-rs blocks",
    )
    p.add_argument("--blocks", required=True, help="IBD blocks TSV (hmmibd-rs output)")
    p.add_argument("--snps", required=True, help="SNP file (VCF or BED)")
    p.add_argument("--snp-format", choices=["vcf", "bed"], default="vcf",
                   help="Format of SNP file (default: vcf)")
    p.add_argument("--output", default="ibd_matrix", help="Output prefix (default: ibd_matrix)")
    p.add_argument("--all-blocks", action="store_true",
                   help="Include blocks with different>0 (default: only different==0)")
    p.add_argument("--min-block-snp", type=int, default=IBD_MIN_BLOCK_SNP,
                   help="Drop IBD segments with fewer than this many SNPs; short, SNP-poor "
                        "segments are commonly spurious (default: %(default)s, 0 disables)")
    p.add_argument("--min-block-kb", type=float, default=IBD_MIN_BLOCK_KB,
                   help="Drop IBD segments shorter than this many kb "
                        "(default: %(default)s, 0 disables)")
    p.add_argument("--sep", default="\t", help="Separator for blocks file (default: tab)")
    return p


def parse_args_build_matrix():
    return get_parser_build_matrix().parse_args()


def build_matrix():
    args = parse_args_build_matrix()

    print(f"\n--- Loading SNPs from {args.snps} ({args.snp_format}) ---")
    panel = SnpPanel.load(args.snps, args.snp_format)
    print(f"Loaded {len(panel):,} SNPs")

    print(f"\n--- Loading IBD blocks from {args.blocks} ---")
    blocks_df = ibd_matrix.read_blocks(args.blocks, sep=args.sep)
    print(f"Loaded {len(blocks_df):,} blocks")

    pair_to_row, pair_labels = ibd_matrix.build_pair_index(blocks_df)
    print(f"Found {len(pair_labels):,} unique pairs")

    print("\n--- Building IBD matrix ---")
    print(f"Keeping IBD segments with >= {args.min_block_snp} SNPs and "
          f">= {args.min_block_kb} kb")
    mat = ibd_matrix.build_matrix(
        blocks_df, panel, pair_to_row,
        only_different_zero=not args.all_blocks,
        min_block_snp=args.min_block_snp, min_block_kb=args.min_block_kb,
    )
    density = mat.nnz / (mat.shape[0] * mat.shape[1]) if mat.shape[0] and mat.shape[1] else 0.0
    print(f"Matrix: {mat.shape[0]:,} pairs x {mat.shape[1]:,} SNPs  "
          f"density={density:.4%}  ({mat.nnz:,} non-zero)")

    print(f"\n--- Saving to {args.output}.* ---")
    ibd_matrix.save_matrix(mat, pair_labels, panel.labels, args.output)

    summary = ibd_matrix.pair_summary(mat, pair_labels)
    summary.to_csv(f"{args.output}.pair_summary.csv", index=False)
    print(f"Per-pair summary -> {args.output}.pair_summary.csv")


if __name__ == "__main__":
    build_matrix()
