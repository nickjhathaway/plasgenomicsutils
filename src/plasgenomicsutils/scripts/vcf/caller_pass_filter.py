#!/usr/bin/env python
"""Keep the records the caller itself passed, counting what its FILTER flags remove."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_caller_pass_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils caller_pass_filter",
        description="Keep records whose FILTER is PASS or '.', dropping everything the "
                    "caller itself flagged -- VQSR tranches, Low_VQSLOD, region classes, "
                    "MissingVQSLOD on contigs VQSR never scored -- and counting the removals "
                    "by flag. Its own step so that verdict stays apart from the metric "
                    "thresholds hard_qc_filter applies, which do not act on the FILTER "
                    "column at all.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--allow", nargs="*", default=[], metavar="FLAG",
                   help="FILTER flags to tolerate: a record whose flags are all listed "
                        "here is kept. E.g. --allow MissingVQSLOD Mitochondrion Apicoplast "
                        "keeps the organelle records a nuclear VQSR model could not score.")
    p.add_argument("--keep-bed", default=None,
                   help="Whitelist BED of regions to keep whatever this filter says "
                        "(0-based half-open). Whitelisted variants still face every other "
                        "filter -- this only exempts them from this one.")
    return p


def parse_args_caller_pass_filter():
    return get_parser_caller_pass_filter().parse_args()


def caller_pass_filter():
    args = parse_args_caller_pass_filter()
    F.caller_pass_filter(args.input, args.output, allow=args.allow, keep_bed=args.keep_bed)
    report_counts(args.input, args.output, "caller_pass_filter")


if __name__ == "__main__":
    caller_pass_filter()
