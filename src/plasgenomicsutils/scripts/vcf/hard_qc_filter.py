#!/usr/bin/env python
"""GATK-style hard filter on INFO metrics; keep PASS."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_hard_qc_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils hard_qc_filter",
        description="Hard filter on QD/MQ/SOR/MQRankSum/ReadPosRankSum, keep PASS.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--qd", type=float, default=20.0, help="Drop QD < this (default: 20)")
    p.add_argument("--mq", type=float, default=55.0, help="Drop MQ < this (default: 55)")
    p.add_argument("--sor", type=float, default=3.0, help="Drop SOR > this (default: 3)")
    p.add_argument("--mqranksum", type=float, default=-5.0, help="Drop MQRankSum < this (default: -5)")
    p.add_argument("--readposranksum", type=float, default=-5.0,
                   help="Drop ReadPosRankSum < this (default: -5)")
    p.add_argument("--fs", type=float, default=None, help="Optional: drop FS > this")
    p.add_argument("--keep-bed", default=None,
                   help="Whitelist BED of regions to keep whatever this filter says "
                        "(0-based half-open, like any BED). Whitelisted variants still "
                        "face every other filter -- this only exempts them from this "
                        "one.")
    return p


def parse_args_hard_qc_filter():
    return get_parser_hard_qc_filter().parse_args()


def hard_qc_filter():
    args = parse_args_hard_qc_filter()
    F.hard_qc_filter(args.input, args.output, qd=args.qd, mq=args.mq, sor=args.sor,
                     mqranksum=args.mqranksum, readposranksum=args.readposranksum, fs=args.fs,
                 keep_bed=args.keep_bed)
    report_counts(args.input, args.output, "hard_qc_filter")


if __name__ == "__main__":
    hard_qc_filter()
