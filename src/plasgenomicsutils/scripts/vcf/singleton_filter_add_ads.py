#!/usr/bin/env python
"""Drop near-private variants and add the FORMAT/ADS summed-depth tag."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_singleton_filter_add_ads() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils singleton_filter_add_ads",
        description="Drop variants ALT in <= --min-samples samples; add FORMAT/ADS.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-samples", type=int, default=1,
                   help="Keep variants with more than this many non-ref, non-missing "
                        "genotype calls (default: 1, i.e. drop singletons)")
    return p


def parse_args_singleton_filter_add_ads():
    return get_parser_singleton_filter_add_ads().parse_args()


def singleton_filter_add_ads():
    args = parse_args_singleton_filter_add_ads()
    F.singleton_add_ads(args.input, args.output, min_samples=args.min_samples)
    report_counts(args.input, args.output, "singleton_filter_add_ads")


if __name__ == "__main__":
    singleton_filter_add_ads()
