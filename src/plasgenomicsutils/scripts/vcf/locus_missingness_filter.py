#!/usr/bin/env python
"""Keep loci with low missingness and high per-sample coverage."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_locus_missingness_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils locus_missingness_filter",
        description="Keep loci with F_MISSING < max AND >= frac of samples at ADS >= min.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--f-missing-max", type=float, default=0.05,
                   help="Maximum fraction of missing genotypes (default: 0.05)")
    p.add_argument("--ads-min", type=int, default=10, help="Per-sample coverage threshold (default: 10)")
    p.add_argument("--sample-frac-min", type=float, default=0.95,
                   help="Minimum fraction of samples at ADS >= --ads-min (default: 0.95)")
    p.add_argument("--keep-bed", default=None,
                   help="Whitelist BED of regions to keep whatever this filter says "
                        "(0-based half-open, like any BED). Whitelisted variants still "
                        "face every other filter -- this only exempts them from this "
                        "one.")
    return p


def parse_args_locus_missingness_filter():
    return get_parser_locus_missingness_filter().parse_args()


def locus_missingness_filter():
    args = parse_args_locus_missingness_filter()
    F.locus_missingness_filter(args.input, args.output, f_missing_max=args.f_missing_max,
                               ads_min=args.ads_min, sample_frac_min=args.sample_frac_min,
                 keep_bed=args.keep_bed)
    report_counts(args.input, args.output, "locus_missingness_filter")


if __name__ == "__main__":
    locus_missingness_filter()
