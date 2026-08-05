#!/usr/bin/env python
"""Adjacency list of the sample pairs sharing IBD over each gene."""

from __future__ import annotations

import argparse

from ...lib.ibd_gene_pairs import gene_ibd_pairs
from ...lib.ibd_matrix import IBD_MIN_BLOCK_KB, IBD_MIN_BLOCK_SNP, read_blocks
from ...utils.small_utils import Utils


def get_parser_gene_pairs() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils ibd_gene_pairs",
        description="Sample pairs that are IBD over each gene, with how much of the gene "
                    "each IBD block covers (one table for all the genes given)",
    )
    p.add_argument("--blocks", required=True,
                   help="hmmibd-rs blocks TSV (sample1, sample2, chr, start, end, different)")
    p.add_argument("--genes", required=True,
                   help="Gene intervals: columns name, chr (or chrom), start, end "
                        "(gene_id optional); several genes come back in one table")
    p.add_argument("--gene", action="append", default=None, metavar="NAME",
                   help="Restrict to this gene name from --genes; repeatable")
    p.add_argument("--within", type=int, default=0,
                   help="Pad each gene by this many bp when deciding overlap; coverage is "
                        "still measured against the gene itself (default: 0)")
    p.add_argument("--complete-only", action="store_true",
                   help="Keep only pairs whose block spans the whole gene")
    p.add_argument("--min-percent-covered", type=float, default=None, metavar="PCT",
                   help="Keep only rows covering at least this percent of the gene")
    p.add_argument("--min-block-snp", type=int, default=IBD_MIN_BLOCK_SNP,
                   help="Drop IBD segments with fewer than this many SNPs; short, SNP-poor "
                        "segments are commonly spurious (default: %(default)s, 0 disables)")
    p.add_argument("--min-block-kb", type=float, default=IBD_MIN_BLOCK_KB,
                   help="Drop IBD segments shorter than this many kb "
                        "(default: %(default)s, 0 disables)")
    p.add_argument("--output", required=True, help="Output TSV(.gz)")
    return p


def parse_args_gene_pairs():
    return get_parser_gene_pairs().parse_args()


def gene_pairs():
    args = parse_args_gene_pairs()

    genes = Utils.read_table(args.genes)
    if "chr" not in genes.columns and "chrom" in genes.columns:
        genes = genes.rename(columns={"chrom": "chr"})
    missing = {"name", "chr", "start", "end"} - set(genes.columns)
    if missing:
        raise SystemExit(f"--genes is missing column(s): {', '.join(sorted(missing))}")
    if args.gene:
        want = {g.lower() for g in args.gene}
        unknown = want - set(genes["name"].astype(str).str.lower())
        if unknown:
            raise SystemExit(f"--gene not found in --genes: {', '.join(sorted(unknown))}")
        genes = genes[genes["name"].astype(str).str.lower().isin(want)]

    print(f"Loading blocks from {args.blocks} ...")
    blocks = read_blocks(args.blocks)
    print(f"  {len(blocks):,} segments; listing IBD pairs over {len(genes):,} gene(s) "
          f"(within={args.within}) ...")
    out = gene_ibd_pairs(blocks, genes, within=args.within,
                         min_block_snp=args.min_block_snp,
                         min_block_kb=args.min_block_kb)

    if args.complete_only:
        out = out[out["coverage"] == "complete"]
    if args.min_percent_covered is not None:
        out = out[out["percent_covered"] >= args.min_percent_covered]

    Utils.write_tsv_gz(out, args.output)
    n_genes = out["gene"].nunique() if len(out) else 0
    n_pairs = len(out.drop_duplicates(["sample1", "sample2"])) if len(out) else 0
    print(f"  -> {args.output}  ({len(out):,} rows, {n_pairs:,} distinct pairs, "
          f"{n_genes} gene(s))")


if __name__ == "__main__":
    gene_pairs()
