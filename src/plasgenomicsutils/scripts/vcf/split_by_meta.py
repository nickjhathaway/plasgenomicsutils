#!/usr/bin/env python
"""Split a callset into one VCF per metadata group (e.g. per country)."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F


def get_parser_split_by_meta() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils split_by_meta",
        description="Split samples into one VCF per metadata group; refill AC/AN/AF/MAF, "
                    "with an optional per-group MAF floor and optional ALT trimming.",
    )
    p.add_argument("--input", required=True, help="Input VCF/BCF (all groups together).")
    p.add_argument("--outdir", required=True,
                   help="Directory for the per-group outputs (created if needed). Files are "
                        "named <group><ext>, ext from --output-type.")
    p.add_argument("--meta", required=True,
                   help="Per-sample metadata TSV (with a header), same shape as maf_filter's.")
    p.add_argument("--group-col", required=True,
                   help="Metadata column with the group label to split on (e.g. country).")
    p.add_argument("--sample-col", default="sample",
                   help="Metadata column with the sample id (default: sample).")
    p.add_argument("--maf-min", type=float, default=None,
                   help="Optional per-group minor-allele-frequency floor: within each output, "
                        "keep only sites with MAF >= this (implies --refill).")
    refill = p.add_mutually_exclusive_group()
    refill.add_argument("--refill", dest="refill", action="store_true", default=True,
                        help="Recompute AC,AN,AF,MAF on each subset (default).")
    refill.add_argument("--no-refill", dest="refill", action="store_false",
                        help="Leave the parent's tags as-is (stale after subsetting).")
    p.add_argument("--trim-alts", action="store_true",
                   help="Drop ALT alleles unused within a group (bcftools --trim-alt-alleles). "
                        "Off by default so the pieces remain losslessly bcftools-merge-able.")
    p.add_argument("--output-type", choices=["b", "z", "v"], default="b",
                   help="bcftools -O letter for the outputs: b BCF (default), z vcf.gz, v VCF.")
    return p


def parse_args_split_by_meta():
    return get_parser_split_by_meta().parse_args()


def split_by_meta():
    args = parse_args_split_by_meta()
    F.split_by_meta(args.input, args.outdir, meta=args.meta, group_col=args.group_col,
                    sample_col=args.sample_col, maf_min=args.maf_min, refill=args.refill,
                    trim_alts=args.trim_alts, output_type=args.output_type)


if __name__ == "__main__":
    split_by_meta()
