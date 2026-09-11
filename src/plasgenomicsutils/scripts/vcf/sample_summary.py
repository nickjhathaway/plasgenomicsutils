#!/usr/bin/env python
"""Per-sample coverage and Fws for a callset as it stands; a table, no filtering."""

from __future__ import annotations

import argparse

from ...lib.callset_summary import (sample_summary_note, sample_summary_table,
                                    write_sample_summary)
from ...lib.reporting import say


def get_parser_sample_summary() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils sample_summary",
        description="One row per sample: the coverage columns sample_coverage_filter "
                    "decides on and the Fws columns fws_filter decides on, side by side, "
                    "with nothing dropped. would_drop_coverage / would_drop_fws say what "
                    "those filters would do at the given thresholds. For a final callset "
                    "those filters never saw. FORMAT/ADS is added on a temporary copy when "
                    "the input lacks it.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True, help="TSV to write")
    p.add_argument("--ads-min", type=int, default=10, help="A locus counts as covered at ADS >= this (default 10)")
    p.add_argument("--frac-min", type=float, default=0.80, help="Coverage bar the would_drop column reports against (default 0.80)")
    p.add_argument("--fws-min", type=float, default=0.95, help="Fws bar for monoclonal (default 0.95)")
    p.add_argument("--estimator", default="regression", choices=("regression", "ratio"))
    p.add_argument("--min-depth", type=int, default=0)
    p.add_argument("--n-bins", type=int, default=10)
    p.add_argument("--min-alt-samples", type=int, default=0)
    p.add_argument("--all-variants", dest="snps_only", action="store_false", help="Score Fws over every record, not SNPs only")
    p.add_argument("--multiallelic", default="collapse")
    p.add_argument("--no-trim", dest="trim", action="store_false")
    return p


def sample_summary():
    args = get_parser_sample_summary().parse_args()
    rows, n_sites = sample_summary_table(
        args.input, ads_min=args.ads_min, frac_min=args.frac_min, fws_min=args.fws_min,
        estimator=args.estimator, min_depth=args.min_depth, n_bins=args.n_bins,
        min_alt_samples=args.min_alt_samples, snps_only=args.snps_only,
        multiallelic=args.multiallelic, trim=args.trim)
    write_sample_summary(rows, args.output)
    say(sample_summary_note(rows, n_sites, frac_min=args.frac_min, fws_min=args.fws_min))
    say(f"  -> {args.output}")


if __name__ == "__main__":
    sample_summary()
