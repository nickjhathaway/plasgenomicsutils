#!/usr/bin/env python
"""Keep only variants inside the core-genome BED."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.assets import resolve_bed
from ...lib.bcftools import report_counts


def get_parser_core_region_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils core_region_filter",
        description="Keep only variants within the core genome (drop subtelomeric/hypervariable).",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--bed", default="builtin:pf3d7_core_regions",
                   help="Core-genome BED, or builtin:pf3d7_core_regions (default).")
    p.add_argument("--keep-bed", default=None,
                   help="Whitelist BED of regions to keep whatever this filter says "
                        "(0-based half-open, like any BED). Whitelisted variants still "
                        "face every other filter -- this only exempts them from this "
                        "region rule.")
    return p


def parse_args_core_region_filter():
    return get_parser_core_region_filter().parse_args()


def core_region_filter():
    args = parse_args_core_region_filter()
    F.core_region_filter(args.input, args.output, bed=resolve_bed(args.bed),
                 keep_bed=resolve_bed(args.keep_bed) if args.keep_bed else None)
    report_counts(args.input, args.output, "core_region_filter")


if __name__ == "__main__":
    core_region_filter()
