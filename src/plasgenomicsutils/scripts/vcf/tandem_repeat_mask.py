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
    return p


def parse_args_tandem_repeat_mask():
    return get_parser_tandem_repeat_mask().parse_args()


def tandem_repeat_mask():
    args = parse_args_tandem_repeat_mask()
    F.tandem_repeat_mask(args.input, args.output, bed=resolve_bed(args.bed))
    report_counts(args.input, args.output, "tandem_repeat_mask")


if __name__ == "__main__":
    tandem_repeat_mask()
