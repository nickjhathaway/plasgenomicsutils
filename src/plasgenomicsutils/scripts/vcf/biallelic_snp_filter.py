#!/usr/bin/env python
"""Keep biallelic SNPs (trimming ALT alleles unused after re-genotyping)."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_biallelic_snp_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils biallelic_snp_filter",
        description="Keep biallelic SNPs; trim ALT alleles absent from all genotypes first.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--no-trim", action="store_true",
                   help="Do not trim unused ALT alleles before the biallelic test "
                        "(by default a site left biallelic once an artifact allele is "
                        "re-genotyped away is kept rather than discarded).")
    return p


def parse_args_biallelic_snp_filter():
    return get_parser_biallelic_snp_filter().parse_args()


def biallelic_snp_filter():
    args = parse_args_biallelic_snp_filter()
    F.biallelic_snp_filter(args.input, args.output, trim=not args.no_trim)
    report_counts(args.input, args.output, "biallelic_snp_filter")


if __name__ == "__main__":
    biallelic_snp_filter()
