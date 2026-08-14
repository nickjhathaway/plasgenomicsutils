#!/usr/bin/env python
"""Remove SNPs overlapping a tandem-repeat BED (artifact-prone regions)."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.assets import resolve_bed
from ...lib.bcftools import report_counts


def get_parser_tandem_repeat_mask() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils tandem_repeat_mask",
        description="Remove variants inside tandem-repeat regions (bedtools intersect -v).",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--bed", default="builtin:pf3d7_tandem_repeats",
                   help="Tandem-repeat BED to exclude, or builtin:pf3d7_tandem_repeats (default).")
    p.add_argument("--keep-bed", default=None,
                   help="Whitelist BED of regions to keep whatever this filter says "
                        "(0-based half-open, like any BED). Whitelisted variants still "
                        "face every other filter -- this only exempts them from this "
                        "region rule.")
    return p


def parse_args_tandem_repeat_mask():
    return get_parser_tandem_repeat_mask().parse_args()


def tandem_repeat_mask():
    args = parse_args_tandem_repeat_mask()
    F.tandem_repeat_mask(args.input, args.output, bed=resolve_bed(args.bed),
                 keep_bed=resolve_bed(args.keep_bed) if args.keep_bed else None)
    report_counts(args.input, args.output, "tandem_repeat_mask")


if __name__ == "__main__":
    tandem_repeat_mask()
