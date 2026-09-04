#!/usr/bin/env python
"""Keep SNPs and/or biallelic sites (trimming ALT alleles unused after re-genotyping)."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts


def get_parser_biallelic_snp_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils biallelic_snp_filter",
        description="Keep biallelic SNPs; trim ALT alleles absent from all genotypes first. "
                    "The two tests are separate: --no-biallelic keeps multiallelic SNPs for "
                    "downstream tools that read them, and --no-snps-only keeps indels.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--no-trim", action="store_true",
                   help="Do not trim unused ALT alleles before the biallelic test "
                        "(by default a site left biallelic once an artifact allele is "
                        "re-genotyped away is kept rather than discarded).")
    p.add_argument("--no-snps-only", dest="snps_only", action="store_false",
                   help="Keep every variant type, not only substitutions. SNP selection is "
                        "by exclusion of the other types, so a mixed SNP+indel record does "
                        "not slip through the way `bcftools view -v snps` lets it.")
    p.add_argument("--no-biallelic", dest="biallelic", action="store_false",
                   help="Keep sites with more than one ALT allele. The rest of the site set "
                        "is unchanged; this only stops multiallelic SNPs being dropped.")
    p.add_argument("--mnp-handling", choices=list(F.MNP_HANDLING), default="split",
                   help="What to do with a multi-base substitution: `split` (default) "
                        "breaks it into its component SNPs with `bcftools norm -a`, which "
                        "also rewrites a SNP that was merely written with padding "
                        "(REF=TTATA ALT=CTATA differs at one base) into its minimal form; "
                        "`remove` drops the record; `keep` leaves it for a downstream tool "
                        "that reads MNPs.")
    return p


def parse_args_biallelic_snp_filter():
    return get_parser_biallelic_snp_filter().parse_args()


def biallelic_snp_filter():
    args = parse_args_biallelic_snp_filter()
    F.biallelic_snp_filter(args.input, args.output, trim=not args.no_trim,
                           snps_only=args.snps_only, biallelic=args.biallelic,
                           mnp_handling=args.mnp_handling)
    report_counts(args.input, args.output, "biallelic_snp_filter")


if __name__ == "__main__":
    biallelic_snp_filter()
