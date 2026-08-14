#!/usr/bin/env python
"""Drop low-coverage samples, then refresh AC/AN/AF and re-drop non-variant sites."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_sample_coverage_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils sample_coverage_filter",
        description="Drop samples covered (ADS >= --ads-min) at < --frac-min of loci.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--ads-min", type=int, default=10, help="Coverage threshold on ADS (default: 10)")
    p.add_argument("--frac-min", type=float, default=0.80,
                   help="Minimum fraction of loci a sample must cover to be kept (default: 0.80)")
    p.add_argument("--dropped-samples", default=None,
                   help="Optional path to write the list of dropped samples")
    p.add_argument("--cov-table", default=None,
                   help="Write the per-sample coverage table the keep/drop decision is read "
                        "from (sample, n_covered, frac_covered, mean_ads, margin, dropped). "
                        "filter_pipeline writes this beside each step by default.")
    return p


def parse_args_sample_coverage_filter():
    return get_parser_sample_coverage_filter().parse_args()


def sample_coverage_filter():
    args = parse_args_sample_coverage_filter()
    dropped = F.sample_coverage_filter(
        args.input, args.output, ads_min=args.ads_min, frac_min=args.frac_min,
        dropped_samples_path=args.dropped_samples, cov_table_path=args.cov_table)
    print(f"  dropped {len(dropped)} low-coverage sample(s)"
          + (f": {', '.join(dropped)}" if dropped else "")
          + (f"\n  coverage table -> {args.cov_table}" if args.cov_table else ""))
    report_counts(args.input, args.output, "sample_coverage_filter")


if __name__ == "__main__":
    sample_coverage_filter()
