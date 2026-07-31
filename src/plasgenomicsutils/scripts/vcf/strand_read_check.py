#!/usr/bin/env python
"""Read-level strand-bias diagnostic for one BAM at one position (+ optional read dump)."""

from __future__ import annotations

import argparse
import os
import sys

from ...lib.pileup import lookup_ref_base, parse_pos
from ...lib.strandbias import (
    PER_READ_COLUMNS, extract_alt_reads, format_read_check_report, per_read_records,
)


def get_parser_strand_read_check() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils strand_read_check",
        description="Tabulate every read over a site (strand, base, BQ, true sequencing "
                    "cycle, soft-clip, SA tag, priming start) and summarize the "
                    "reverse-strand ALT reads where SSE artifacts live.")
    p.add_argument("--bam", required=True)
    p.add_argument("--pos", required=True, help="CHROM:POS (1-based)")
    p.add_argument("--ref", help="Reference FASTA (to look up the ref base and, with "
                                 "--extract-reads, pad the aligned FASTA)")
    p.add_argument("--ref-base", help="Reference base at POS (overrides --ref lookup)")
    p.add_argument("--alt-base", help="Focus one ALT base (default: any non-ref)")
    p.add_argument("--out", help="Output prefix (default: BAM basename + position)")
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument("--extract-reads", action="store_true",
                   help="Also write <prefix>.alt_reads.fastq and .refaln.fasta for viewing")
    p.add_argument("--window", type=int, default=40,
                   help="+/- bp for the aligned FASTA when --extract-reads (default: 40)")
    p.add_argument("--strand", choices=["rev", "fwd", "both"], default="rev",
                   help="Which strand's ALT reads to extract (default: rev = the artifact)")
    p.add_argument("--include-ref", type=int, default=0,
                   help="Include up to N same-strand reference reads for contrast when extracting")
    return p


def parse_args_strand_read_check():
    return get_parser_strand_read_check().parse_args()


def strand_read_check():
    import pysam

    args = parse_args_strand_read_check()
    chrom, pos1 = parse_pos(args.pos)
    pos0 = pos1 - 1

    ref_base = args.ref_base
    if ref_base is None and args.ref:
        ref_base = lookup_ref_base(args.ref, chrom, pos0)
    if ref_base is None:
        sys.exit("Provide --ref-base or --ref so the reference allele is known.")
    ref_base = ref_base.upper()

    prefix = args.out or f"{os.path.basename(args.bam).replace('.bam', '')}_{chrom}_{pos1}"
    bam = pysam.AlignmentFile(args.bam, "rb")
    try:
        rows = per_read_records(bam, chrom, pos0, ref_base, min_mapq=args.min_mapq)

        tsv = prefix + ".per_read.tsv"
        with open(tsv, "w") as fh:
            fh.write("\t".join(PER_READ_COLUMNS) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in PER_READ_COLUMNS) + "\n")

        print(format_read_check_report(rows, os.path.basename(args.bam), chrom, pos1,
                                       ref_base, args.alt_base))
        print(f"\n  wrote {tsv}")

        if args.extract_reads:
            if not args.ref:
                sys.exit("--extract-reads needs --ref (the FASTA) to anchor the aligned output.")
            n_alt, n_ref, fq, fa = extract_alt_reads(
                bam, args.ref, chrom, pos1, prefix, alt_base=args.alt_base,
                window=args.window, strand=args.strand, include_ref=args.include_ref,
                min_mapq=args.min_mapq)
            if n_alt == 0:
                print("  no ALT reads to extract with the given --strand/--alt-base filters.")
            else:
                print(f"  extracted {n_alt} ALT reads ({n_ref} ref for contrast):")
                print(f"    {fq}")
                print(f"    {fa}   (POS {pos1} is column {args.window}; reference is first record)")
    finally:
        bam.close()


if __name__ == "__main__":
    strand_read_check()
