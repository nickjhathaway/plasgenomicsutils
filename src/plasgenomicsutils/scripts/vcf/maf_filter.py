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
    p.add_argument("--maf-max", type=float, default=0.99, help="Upper AF bound (bcftools -Q, default: 0.99)")
    return p


def parse_args_maf_filter():
    return get_parser_maf_filter().parse_args()


def maf_filter():
    args = parse_args_maf_filter()
    F.maf_filter(args.input, args.output, maf_min=args.maf_min, maf_max=args.maf_max)
    report_counts(args.input, args.output, "maf_filter")


if __name__ == "__main__":
    maf_filter()
