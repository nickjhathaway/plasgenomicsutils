#!/usr/bin/env python
"""Scan a biallelic VCF/BCF with ADF/ADR for strand-bias (SSE) fake-het artifacts."""

from __future__ import annotations

import argparse
import sys

from ...lib.strandbias import scan_vcf_strand_bias

_TSV_COLUMNS = ["chrom", "pos", "ref", "alt", "sample", "ref_fwd", "alt_fwd", "ref_rev",
                "alt_rev", "vaf_fwd", "vaf_rev", "vaf", "ratio", "sb_phred", "drop",
                "reasons"]


def get_parser_strand_bias_scan() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils strand_bias_scan",
        description="Per-(site, sample) strand-bias verdicts from FORMAT/ADF+ADR. "
                    "Input must be biallelic (bcftools norm -m-) with ADF/ADR present "
                    "(bcftools mpileup -a FORMAT/ADF,FORMAT/ADR).")
    p.add_argument("--input-vcf", required=True, help="Biallelic VCF/BCF with FORMAT/ADF,ADR")
    p.add_argument("--out-tsv", default="-",
                   help="Per-(site, sample) verdict table ('-' = STDOUT, default)")
    p.add_argument("--out-bed",
                   help="Also write a BED of positions dropped in >= --min-drop-samples samples")
    p.add_argument("--min-drop-samples", type=int, default=1,
                   help="Blacklist a position in the BED once this many samples drop it (default: 1)")
    p.add_argument("--min-vaf", type=float, default=0.05,
                   help="Only judge samples whose pooled alt VAF is >= this (default: 0.05)")
    p.add_argument("--min-minor-depth", type=int, default=20,
                   help="Minor strand must have >= this depth for the strand-restricted rule (default: 20)")
    p.add_argument("--ratio-hard", type=float, default=0.15,
                   help="Drop if minor/major strand VAF ratio < this (default: 0.15)")
    p.add_argument("--ratio-soft", type=float, default=0.30,
                   help="Softer ratio paired with low alt BQ (read-level only; default: 0.30)")
    p.add_argument("--sb-hard", type=float, default=60.0,
                   help="Drop if Fisher strand-bias Phred > this (default: 60)")
    return p


def parse_args_strand_bias_scan():
    return get_parser_strand_bias_scan().parse_args()


def strand_bias_scan():
    args = parse_args_strand_bias_scan()
    dropped_positions: dict[tuple, list] = {}

    out = sys.stdout if args.out_tsv == "-" else open(args.out_tsv, "w")
    n_rows = n_drop = 0
    try:
        out.write("\t".join(_TSV_COLUMNS) + "\n")
        for r in scan_vcf_strand_bias(
                args.input_vcf, min_vaf=args.min_vaf, min_minor_depth=args.min_minor_depth,
                ratio_hard=args.ratio_hard, ratio_soft=args.ratio_soft, sb_hard=args.sb_hard):
            n_rows += 1
            out.write("\t".join(_fmt(r[c]) for c in _TSV_COLUMNS) + "\n")
            if r["drop"]:
                n_drop += 1
                key = (r["chrom"], r["pos"], r["ref"], r["alt"])
                dropped_positions.setdefault(key, []).append(r["sample"])
    finally:
        if out is not sys.stdout:
            out.close()

    n_bed = 0
    if args.out_bed:
        with open(args.out_bed, "w") as bed:
            for (chrom, pos, ref, alt), samples in sorted(dropped_positions.items()):
                if len(samples) >= args.min_drop_samples:
                    n_bed += 1
                    bed.write(f"{chrom}\t{pos - 1}\t{pos}\t{ref}>{alt}\t{len(samples)}\n")

    print(f"Done. {n_rows} candidate (site, sample) pairs judged, {n_drop} flagged to drop"
          + (f"; {n_bed} positions written to {args.out_bed}" if args.out_bed else ""),
          file=sys.stderr)


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


if __name__ == "__main__":
    strand_bias_scan()
