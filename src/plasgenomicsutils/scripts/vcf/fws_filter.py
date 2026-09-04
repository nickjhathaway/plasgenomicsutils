#!/usr/bin/env python
"""Keep only monoclonal samples (Fws >= a threshold), then refresh AC/AN/AF."""

from __future__ import annotations

import argparse

from ...lib.bcftools import report_counts
from ...lib.fws import fws_filter as _fws_filter


def get_parser_fws_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils fws_filter",
        description="Drop samples with Fws below --fws-min, keeping the monoclonal "
                    "infections. Samples are dropped and no variants are: removing samples "
                    "changes every allele frequency, so re-run maf_filter / "
                    "locus_missingness_filter afterwards if a site set has to match the "
                    "cohort that remains.",
        epilog="Fws is measured against the cohort's own allele frequencies, so run this on "
               "a callset the rest of the chain has already filtered and re-genotyped --\n"
               "unfiltered calls read as within-host diversity and push every sample down.\n"
               "Use `calculate_fws` when the scores are what you want, not a filtered "
               "callset.\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--fws-min", type=float, default=0.95,
                   help="Keep samples with Fws >= this (default: 0.95, the same bar "
                        "`calculate_fws --monoclonal-threshold` reports on)")
    p.add_argument("--estimator", choices=["regression", "ratio"], default="regression",
                   help="Fws estimator (default: regression, matching moimix::getFws)")
    p.add_argument("--min-depth", type=int, default=0,
                   help="Ignore per-sample sites below this total depth (default: 0)")
    p.add_argument("--n-bins", type=int, default=10, help="Number of MAF bins (default: 10)")
    p.add_argument("--min-alt-samples", type=int, default=0,
                   help="Ignore sites with fewer than this many alt-carrying samples")
    p.add_argument("--no-snps-only", dest="snps_only", action="store_false",
                   help="Score on every biallelic record, not only single-base REF/ALT")
    p.add_argument("--exclude-call-regions", default=None,
                   help="BED of regions to exclude from the Fws calculation")
    p.add_argument("--dropped-samples", default=None,
                   help="Optional path to write the list of dropped samples")
    p.add_argument("--fws-table", default=None,
                   help="Write the per-sample scores the keep/drop decision is read from "
                        "(sample, fws, n_sites, monoclonal, dropped). filter_pipeline "
                        "writes this beside the step by default.")
    return p


def parse_args_fws_filter():
    return get_parser_fws_filter().parse_args()


def fws_filter():
    args = parse_args_fws_filter()
    dropped = _fws_filter(
        args.input, args.output, fws_min=args.fws_min, estimator=args.estimator,
        min_depth=args.min_depth, n_bins=args.n_bins,
        min_alt_samples=args.min_alt_samples, snps_only=args.snps_only,
        exclude_call_regions=args.exclude_call_regions,
        dropped_samples_path=args.dropped_samples, fws_table_path=args.fws_table)
    print(f"  dropped {len(dropped)} sample(s) below Fws {args.fws_min:g}"
          + (f": {', '.join(dropped)}" if dropped else "")
          + (f"\n  Fws table -> {args.fws_table}" if args.fws_table else ""))
    report_counts(args.input, args.output, "fws_filter")


if __name__ == "__main__":
    fws_filter()
