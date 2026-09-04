#!/usr/bin/env python
"""Drop records with no ALT allele, counting them as their own step."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_no_alt_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils no_alt_filter",
        description="Remove records with no ALT allele -- positions that turned out "
                    "non-variant, which calling a region list without -v reports. Its own "
                    "step so the count is explicit: how many positions had nothing to "
                    "call, kept apart from how many real variants failed QC.",
        epilog="Worth knowing: bcftools computes FS, RPBZ, MQBZ and the rest whether or "
               "not an ALT was called, so hard_qc_filter does remove non-variant records "
               "on its own. Running this first is what makes that visible.\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--keep-no-alts", action="store_true",
                   help="Keep the non-variant records instead of dropping them, for a "
                        "fill-in workflow where 'this sample is reference here' is the "
                        "answer being sought. They are counted either way.")
    p.add_argument("--trim-alt-alleles", dest="trim", action="store_true",
                   help="First drop ALT alleles no genotype carries. A callset subset from "
                        "a larger cohort keeps every ALT the full cohort had, now carried "
                        "by nobody -- `bcftools view -S` removes samples, not alleles -- so "
                        "without this the counts describe a cohort that is not the one in "
                        "the file. Reported separately from the no-ALT count, since "
                        "trimming is what leaves some records with no ALT at all.")
    p.add_argument("--keep-bed", default=None,
                   help="Whitelist BED of regions to keep whatever this filter says "
                        "(0-based half-open). Whitelisted variants still face every other "
                        "filter -- this only exempts them from this one.")
    return p


def parse_args_no_alt_filter():
    return get_parser_no_alt_filter().parse_args()


def no_alt_filter():
    args = parse_args_no_alt_filter()
    F.no_alt_filter(args.input, args.output, keep=args.keep_no_alts, trim=args.trim,
                    keep_bed=args.keep_bed)
    report_counts(args.input, args.output, "no_alt_filter")


if __name__ == "__main__":
    no_alt_filter()
