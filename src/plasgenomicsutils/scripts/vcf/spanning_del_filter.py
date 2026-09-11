#!/usr/bin/env python
"""Recode `*` calls as missing and drop the allele: `*` is missingness, not an allele."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts
from ...lib.reporting import say
from ...lib.spanning_del import spanning_del_note


def get_parser_spanning_del_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils spanning_del_filter",
        description="A `*` in ALT says a deletion called somewhere else covers this "
                    "position in some samples. This sets the calls that name it to missing "
                    "and then drops the allele, leaving the site and its real alleles in "
                    "place -- so the remaining variants are those of the non-deleted "
                    "strains. Opt-in: see the note below on what it costs.",
        epilog="OFF BY DEFAULT in the shipped pipeline config, and deliberately so. A `*` "
               "is a confident observation that the sequence is absent, not a failure to "
               "call: a site with 20 deleted samples and 5 carrying a variant comes out of "
               "this reading as though 20 samples could not be genotyped there. In P. "
               "falciparum that matters more than usual, since across the dimorphic regions "
               "a deletion is often the other haplotype rather than a dropout, and which "
               "samples carry it is a result.\n\n"
               "Turn it on when the question is about the variants of the NON-DELETED "
               "strains and the deletion itself is not the subject. Then the `*` is in the "
               "way: left alone it is scored as a third allele, and a record like `A > *,T` "
               "-- an ordinary SNP that merely sits under someone else's deletion -- looks "
               "multiallelic, so `-M2` deletes it. Run this BEFORE biallelic_snp_filter; "
               "afterwards the records it would rescue are already gone. On the shipped Pf7 "
               "fixture that is 312 records rescued, and the genuine multiallelic count "
               "falls from 544 to 232.\n\n"
               "The cost is reported per run: at the 369 fixture records where `*` sits "
               "beside real alleles, a mean 17.2%% of called samples carry it. A call is "
               "nulled slot by slot, so `*/T` becomes `./T` and keeps the T it does carry; "
               "only a call naming nothing but `*` goes fully missing. Run "
               "locus_missingness_filter after this, not before.\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--no-trim", dest="trim", action="store_false",
                   help="Recode the genotypes but leave `*` in the ALT column. The record "
                        "still looks multiallelic until something trims it; useful for "
                        "checking what the recode alone did. Note the trim also removes "
                        "any other ALT no genotype carries, as biallelic_snp_filter does.")
    return p


def parse_args_spanning_del_filter():
    return get_parser_spanning_del_filter().parse_args()


def spanning_del_filter():
    args = parse_args_spanning_del_filter()
    st = F.spanning_del_filter(args.input, args.output, trim=args.trim)
    say(spanning_del_note(st))
    report_counts(args.input, args.output, "spanning_del_filter")


if __name__ == "__main__":
    spanning_del_filter()
