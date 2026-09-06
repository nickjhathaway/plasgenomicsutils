#!/usr/bin/env python
"""Convert a VCF/BCF to BED, 0-based, on stdout unless --out is given."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F


def get_parser_vcf_to_bed() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils vcf_to_bed",
        description="Write the records of a VCF/BCF as BED: chrom, 0-based start, end, and "
                    "a chrom:pos0 name. VCF POS is 1-based and BED is not, so this is the "
                    "conversion rather than a reformat; the interval spans the REF allele, "
                    "one base for a SNP and len(REF) for an indel.",
        epilog="Examples:\n\n"
               "    plasgenomicsutils vcf_to_bed --input calls.bcf > calls.bed\n"
               "    plasgenomicsutils vcf_to_bed --input calls.bcf --snps-only --out snps.bed\n"
               "    plasgenomicsutils vcf_to_bed --input calls.bcf --no-name | \\\n"
               "        bcftools view -R - other.bcf\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True, help="VCF/BCF (any format bcftools reads)")
    p.add_argument("--out", default=None,
                   help="Write here instead of stdout")
    p.add_argument("--snps-only", action="store_true",
                   help="Keep single-base substitutions only, dropping indels and the rest")
    p.add_argument("--no-name", dest="name_column", action="store_false",
                   help="Emit a bare 3-column BED, without the chrom:pos0 name column")
    return p


def parse_args_vcf_to_bed():
    return get_parser_vcf_to_bed().parse_args()


def vcf_to_bed():
    args = parse_args_vcf_to_bed()
    F.vcf_to_bed(args.input, args.out, snps_only=args.snps_only,
                 name_column=args.name_column)


if __name__ == "__main__":
    vcf_to_bed()
