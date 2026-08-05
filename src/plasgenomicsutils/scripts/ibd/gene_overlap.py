#!/usr/bin/env python
"""Per-gene IBD-block overlap between sample groups (for the R gene triangles)."""

from __future__ import annotations

import argparse

from ...lib.ibd_gene_overlap import gene_block_overlap
from ...lib.ibd_matrix import read_blocks
from ...utils.small_utils import Utils


def get_parser_gene_overlap() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils ibd_gene_overlap",
        description="Fraction of sample pairs whose IBD block overlaps each gene, per group pair",
    )
    p.add_argument("--blocks", required=True,
                   help="hmmibd-rs blocks TSV (sample1, sample2, chr, start, end, different)")
    p.add_argument("--genes", required=True,
                   help="Gene intervals: columns name, chr (or chrom), start, end")
    p.add_argument("--meta", required=True, help="Sample metadata: sample + <group-col>")
    p.add_argument("--group-col", default="group",
                   help="Column in --meta to group samples by (default: group)")
    p.add_argument("--within", type=int, default=0,
                   help="Pad each gene interval by this many bp on both sides (default: 0)")
    p.add_argument("--output", required=True, help="Output TSV(.gz)")
    return p


def parse_args_gene_overlap():
    return get_parser_gene_overlap().parse_args()


def gene_overlap():
    args = parse_args_gene_overlap()

    genes = Utils.read_table(args.genes)
    if "chr" not in genes.columns and "chrom" in genes.columns:
        genes = genes.rename(columns={"chrom": "chr"})
    missing = {"name", "chr", "start", "end"} - set(genes.columns)
    if missing:
        raise SystemExit(f"--genes is missing column(s): {', '.join(sorted(missing))}")

    meta = Utils.read_table(args.meta)
    if args.group_col not in meta.columns:
        raise SystemExit(f"--meta has no column '{args.group_col}'")
    sample_to_group = (
        meta.dropna(subset=[args.group_col])
            .set_index("sample")[args.group_col]
            .astype(str)
            .to_dict()
    )

    print(f"Loading blocks from {args.blocks} ...")
    blocks = read_blocks(args.blocks)
    print(f"  {len(blocks):,} segments; computing overlap for {len(genes):,} genes "
          f"(group_col='{args.group_col}', within={args.within}) ...")
    out = gene_block_overlap(blocks, genes, sample_to_group, within=args.within)

    Utils.write_tsv_gz(out, args.output)
    n_genes = out["gene"].nunique() if len(out) else 0
    print(f"  -> {args.output}  ({n_genes} genes x {out['group_a'].nunique()} groups)")


if __name__ == "__main__":
    gene_overlap()
