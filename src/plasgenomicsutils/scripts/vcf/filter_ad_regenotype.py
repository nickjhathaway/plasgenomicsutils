#!/usr/bin/env python
"""Clean within-sample AD artifacts by depth/frequency, then re-genotype."""

from __future__ import annotations

import argparse

from ...lib.regenotype import filter_ad_regenotype as _run


def get_parser_filter_ad_regenotype() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils filter_ad_regenotype",
        description="Filter AD by read depth and within-sample frequency, then re-genotype.",
    )
    p.add_argument("--input-vcf", required=True, help="Input VCF/BCF (may be gzipped)")
    p.add_argument("--output-vcf", required=True, help="Output VCF/BCF")
    p.add_argument("--min-reads", type=int, default=2,
                   help="Zero AD alleles with < this many reads (default: 2)")
    p.add_argument("--min-freq", type=float, default=0.01,
                   help="Zero AD alleles with within-sample freq (AD/ADS) < this (default: 0.01)")
    p.add_argument("--het-min-af", type=float, default=0.2,
                   help="Minimum minor-allele frequency to call a heterozygote (default: 0.2)")
    p.add_argument("--restrict-to-called-alleles", action="store_true",
                   help="Only narrow each sample's existing genotype instead of re-deriving "
                        "it from AD. A call can lose support and go missing but never gain a "
                        "new allele, preserving the upstream caller's genotypes where they "
                        "disagree with raw read counts.")
    p.add_argument("--ploidy", type=int, choices=(1, 2), default=None,
                   help="Output genotype ploidy. Default keeps the conventional diploid "
                        "coding (0/1 = mixed infection). If set, it is validated against the "
                        "input ploidy per record: greater than the input errors (cannot "
                        "promote), less than warns and trims genotype-linked fields (PL/GL) "
                        "to match. 1 = haploid (single best allele), 2 = diploid.")
    p.add_argument("--no-add-ads", action="store_true",
                   help="Do not derive FORMAT/ADS from AD when the callset lacks it. "
                        "Without ADS every record is passed through unfiltered, which "
                        "looks exactly like a filter that had nothing to do -- so by "
                        "default it is derived (the sum of AD, as singleton_filter_add_ads "
                        "computes it).")
    p.add_argument("--singletons", choices=("none", "real", "star", "all"), default="star",
                   help="After re-genotyping, blank the calls naming an allele carried by "
                        "<= --singleton-min-samples samples and trim the allele off: "
                        "`star` (the default) the spanning deletion `*` only, `real` the "
                        "bases and indels only, `all` both, `none` to skip it. "
                        "Re-genotyping moves the carrier counts the earlier singleton step "
                        "judged, so an allele it saw with two carriers can have one here -- "
                        "and a surviving singleton `*` costs the whole record at "
                        "biallelic_snp_filter --snps-only. A `*` with more carriers is left "
                        "for spanning_del_filter.")
    p.add_argument("--singleton-min-samples", type=int, default=1,
                   help="Carrier count at or below which --singletons blanks an allele "
                        "(default: 1, i.e. private to one sample)")
    return p


def parse_args_filter_ad_regenotype():
    return get_parser_filter_ad_regenotype().parse_args()


def filter_ad_regenotype():
    args = parse_args_filter_ad_regenotype()
    _run(args.input_vcf, args.output_vcf,
         min_reads=args.min_reads, min_freq=args.min_freq, het_min_af=args.het_min_af,
         restrict_to_called=args.restrict_to_called_alleles, ploidy=args.ploidy,
         add_ads=not args.no_add_ads, singletons=args.singletons,
         singleton_min_samples=args.singleton_min_samples)
    print("Done.")


if __name__ == "__main__":
    filter_ad_regenotype()
