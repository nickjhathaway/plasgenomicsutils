#!/usr/bin/env python
"""Rewrite GATK's balanced indel pairs as the SNPs they encode."""

from __future__ import annotations

import argparse

from ...lib.bcftools import report_counts
from ...lib.shifted_indels import resolve_shifted_indels


def get_parser_resolve_shifted_indels() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils resolve_shifted_indels",
        description="HaplotypeCaller often writes a multi-base substitution as an insertion "
                    "and a deletion a few bases apart -- its haplotype-to-reference alignment "
                    "scores two gaps above three mismatches -- so a coding change like pfcrt "
                    "CVIET arrives as two indel records plus a SNP, and an SNP-only chain "
                    "drops two thirds of it. The records share a PID (GATK's physical "
                    "phasing). This step groups records by PID, rebuilds each sample's "
                    "haplotypes over the cluster, and where every haplotype is reference-"
                    "length replaces the cluster with one SNP record per changed base. "
                    "Clusters with a net length change are real indels, and clusters that "
                    "would cost a sample its genotype for lack of phase, are left alone.",
        epilog="Derived records carry AD/DP and INFO from the nearest source record and "
               "INFO/SHIFTED naming the sources; PL/GQ/PGT/PID/PS/SB are dropped; AC/AN/AF "
               "are refilled.\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--reference", default=None,
                   help="Reference FASTA (indexed). Default: the header's ##reference line.")
    p.add_argument("--max-unresolved", type=int, default=1,
                   help="How many samples the rewrite may cost. A sample heterozygous at two "
                        "or more of a cluster's records with no PGT to phase them cannot be "
                        "laid out, and a cluster that would cost more than this many such "
                        "samples is left as it was. The trade is very uneven: on a "
                        "249-sample cohort the default of 1 gives 1,953 SNPs for 60 "
                        "genotypes, while the worst clusters cost over 25 samples each for "
                        "a handful of SNPs. Use 0 to never write a missing genotype.")
    p.add_argument("--max-span", type=int, default=50,
                   help="Leave clusters whose reference footprint exceeds this many bp "
                        "alone (default: 50)")
    return p


def resolve_shifted_indels_cmd():
    args = get_parser_resolve_shifted_indels().parse_args()
    resolve_shifted_indels(args.input, args.output, reference=args.reference,
                           max_span=args.max_span, max_unresolved=args.max_unresolved)
    report_counts(args.input, args.output, "resolve_shifted_indels")


if __name__ == "__main__":
    resolve_shifted_indels_cmd()
