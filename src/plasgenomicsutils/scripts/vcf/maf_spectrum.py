#!/usr/bin/env python
"""Where a callset's variation sits on the frequency spectrum, overall and per group."""

from __future__ import annotations

import argparse

from ...lib.callset_summary import (MAF_MARKS, maf_spectrum_note, maf_spectrum_table,
                                    write_maf_spectrum)
from ...lib.reporting import say


def get_parser_maf_spectrum() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils maf_spectrum",
        description="How many records sit at or above each minor-allele-frequency mark, so "
                    "the cost of a 1%%, 2%%, 5%% or 10%% floor is a table rather than a rerun. "
                    "Run it BEFORE the frequency filter: afterwards the sub-floor records "
                    "are gone. With --meta and --group-col the spectrum is also computed "
                    "per group, since a grouped floor keeps a record when any one group "
                    "clears the bar.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True, help="TSV to write")
    p.add_argument("--marks", type=float, nargs="+", default=list(MAF_MARKS), metavar="MAF",
                   help=f"Frequency marks to report (default: {' '.join(str(m) for m in MAF_MARKS)})")
    p.add_argument("--meta", default=None, help="Per-sample metadata TSV, for a per-group spectrum")
    p.add_argument("--group-col", default=None, help="Grouping column in --meta")
    p.add_argument("--sample-col", default="sample", help="Sample-name column in --meta")
    return p


def maf_spectrum():
    args = get_parser_maf_spectrum().parse_args()
    rows = maf_spectrum_table(args.input, marks=tuple(args.marks), meta=args.meta,
                              group_col=args.group_col, sample_col=args.sample_col)
    write_maf_spectrum(rows, args.output)
    say(maf_spectrum_note(rows))
    say(f"  -> {args.output}")


if __name__ == "__main__":
    maf_spectrum()
