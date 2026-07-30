#!/usr/bin/env python
"""Compute the Fws within-host diversity statistic from a VCF/BCF or an AD table."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ...lib.fws import compute_fws, load_exclude_regions, read_ad_table, read_ad_vcf


def get_parser_calculate_fws() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils calculate_fws",
        description="Per-sample Fws (Manske 2012 / moimix::getFws) from per-sample AD, "
                    "read from a VCF/BCF or a bcftools-query AD table.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--input-vcf", help="VCF/BCF with FORMAT/AD (biallelic SNPs used; "
                                         "split multiallelics with `bcftools norm -m-` first)")
    src.add_argument("--ad-table", help="bcftools query TSV: CHROM POS REF ALT then one AD "
                                        "('ref,alt') per sample")
    p.add_argument("--samples", help="For --ad-table: comma-separated sample ids in column "
                                     "order, or a file with one per line")
    p.add_argument("--estimator", choices=["regression", "ratio"], default="regression",
                   help="regression = moimix::getFws (default); ratio = wgs_cnv_workflow legacy")
    p.add_argument("--min-depth", type=int, default=0,
                   help="Drop per-sample sites below this read depth (default: 0; CNV gate uses 10)")
    p.add_argument("--n-bins", type=int, default=10, help="Number of MAF bins (default: 10)")
    p.add_argument("--min-alt-samples", type=int, default=0,
                   help="Keep only sites with the alt seen in >= this many samples "
                        "(default: 0; CNV gate uses 2)")
    p.add_argument("--snps-only", action="store_true",
                   help="Keep only single-base REF/ALT sites (default: all biallelic, "
                        "matching moimix)")
    p.add_argument("--monoclonal-threshold", type=float, default=0.95,
                   help="Fws >= this is reported monoclonal (default: 0.95)")
    p.add_argument("--population-name",
                   help="If set, add a population_name column with this value (eases merging)")
    p.add_argument("--exclude-call-regions",
                   help="TSV (chrom, call_start, call_end) of CNV windows to exclude from Fws")
    p.add_argument("--out", default="-", help="Output TSV ('-' = STDOUT, default)")
    return p


def parse_args_calculate_fws():
    return get_parser_calculate_fws().parse_args()


def _load_samples(spec):
    sp = Path(spec)
    if sp.exists():
        return [x.strip() for x in sp.read_text().splitlines() if x.strip()]
    return [x for x in spec.split(",") if x]


def calculate_fws():
    args = parse_args_calculate_fws()
    exclude = load_exclude_regions(args.exclude_call_regions)

    if args.input_vcf:
        samples, ref, alt = read_ad_vcf(args.input_vcf, exclude, snps_only=args.snps_only)
    else:
        if not args.samples:
            sys.exit("--ad-table requires --samples (the column order)")
        samples = _load_samples(args.samples)
        ref, alt = read_ad_table(args.ad_table, samples, exclude, snps_only=args.snps_only)

    if ref.size == 0:
        sys.exit("no usable biallelic SNP sites found in the input")

    fws, n_info = compute_fws(ref, alt, estimator=args.estimator, min_depth=args.min_depth,
                              n_bins=args.n_bins, min_alt_samples=args.min_alt_samples)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        header = ["sample", "fws", "n_sites", "monoclonal"]
        if args.population_name is not None:
            header.append("population_name")
        out.write("\t".join(header) + "\n")
        for s, f, n in zip(samples, fws, n_info):
            mono = "" if not np.isfinite(f) else str(bool(f >= args.monoclonal_threshold))
            fval = "" if not np.isfinite(f) else f"{f:.6f}"
            row = [s, fval, str(int(n)), mono]
            if args.population_name is not None:
                row.append(args.population_name)
            out.write("\t".join(row) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()

    n_ok = int(np.isfinite(fws).sum())
    n_mono = int(np.nansum(fws >= args.monoclonal_threshold))
    print(f"fws: {n_ok} samples scored over {ref.shape[0]} sites ({args.estimator}); "
          f"{n_mono} monoclonal (Fws>={args.monoclonal_threshold})", file=sys.stderr)


if __name__ == "__main__":
    calculate_fws()
