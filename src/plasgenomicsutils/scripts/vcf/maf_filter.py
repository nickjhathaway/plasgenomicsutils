#!/usr/bin/env python
"""Drop rare and near-fixed alleles by a minor-allele-frequency window."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_maf_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils maf_filter",
        description="Keep variants with allele frequency in [--maf-min, --maf-max].",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--maf-min", type=float, default=0.01, help="Lower AF bound (bcftools -q, default: 0.01)")
    p.add_argument("--maf-max", type=float, default=None,
                   help="Upper AF bound (bcftools -Q). Defaults to 1 - maf_min (a symmetric window).")
    p.add_argument("--meta", default=None,
                   help="Per-sample metadata TSV (with a header) for per-group MAF. A site is "
                        "kept if its minor-allele frequency is >= --maf-min in ANY group, judged "
                        "on the combined VCF so all genotypes are preserved.")
    p.add_argument("--group-col", default=None,
                   help="Metadata column with the group label (e.g. country); enables grouped MAF.")
    p.add_argument("--sample-col", default="sample",
                   help="Metadata column with the sample id (default: sample).")
    return p


def parse_args_maf_filter():
    return get_parser_maf_filter().parse_args()


def maf_filter():
    args = parse_args_maf_filter()
    F.maf_filter(args.input, args.output, maf_min=args.maf_min, maf_max=args.maf_max,
                 meta=args.meta, group_col=args.group_col, sample_col=args.sample_col)
    report_counts(args.input, args.output, "maf_filter")


if __name__ == "__main__":
    maf_filter()
