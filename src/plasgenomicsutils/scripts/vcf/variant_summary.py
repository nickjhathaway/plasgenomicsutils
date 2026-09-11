#!/usr/bin/env python
"""Records by class and by ALT-allele count, as counts and fractions; a table, no filtering."""

from __future__ import annotations

import argparse

from ...lib.callset_summary import (variant_summary_note, variant_summary_table,
                                    write_variant_summary)
from ...lib.reporting import say


def get_parser_variant_summary() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils variant_summary",
        description="Every record classed (snps, indels, mnps, mixed, spanning_del, ...) "
                    "and, within each class, broken down by how many ALT alleles it "
                    "carries -- biallelic, triallelic, ... -- as counts, as a fraction of "
                    "the whole callset and as a fraction of the class.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True, help="TSV to write")
    return p


def variant_summary():
    args = get_parser_variant_summary().parse_args()
    rows = variant_summary_table(args.input)
    write_variant_summary(rows, args.output)
    say(variant_summary_note(rows))
    say(f"  -> {args.output}")


if __name__ == "__main__":
    variant_summary()
