#!/usr/bin/env python
"""Remove variants in paralogous / multigene-family genes (mismapping-prone)."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.assets import resolve_bed
from ...lib.bcftools import report_counts


def get_parser_paralog_mask() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils paralog_mask",
        description="Remove variants overlapping paralogous/multigene-family genes.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--bed", default="builtin:pf3d7_paralog_genes",
                   help="Paralog-gene BED, or builtin:pf3d7_paralog_genes (default).")
    p.add_argument("--keep-bed", default=None,
                   help="Whitelist BED of regions to keep whatever this filter says "
                        "(0-based half-open, like any BED). Whitelisted variants still "
                        "face every other filter -- this only exempts them from this "
                        "region rule.")
    return p


def parse_args_paralog_mask():
    return get_parser_paralog_mask().parse_args()


def paralog_mask():
    args = parse_args_paralog_mask()
    F.paralog_mask(args.input, args.output, bed=resolve_bed(args.bed),
                 keep_bed=resolve_bed(args.keep_bed) if args.keep_bed else None)
    report_counts(args.input, args.output, "paralog_mask")


if __name__ == "__main__":
    paralog_mask()
