#!/usr/bin/env python
"""Per-chromosome distances between consecutive variants, and variant density."""

from __future__ import annotations

import argparse

from ...lib.reference import DEFAULT_REFERENCE, available_references, get_reference
from ...lib.variant_spacing import variant_spacing as _spacing
from ...lib.variant_spacing import write_table


def get_parser_variant_spacing() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils variant_spacing",
        description="How far apart the variants are, per chromosome: the quartiles of the "
                    "gap between consecutive variants, plus density per cM. Density alone "
                    "cannot tell an evenly spaced chromosome from dense clusters separated "
                    "by deserts; the spread of the gaps can.",
        epilog="Examples:\n\n"
               "    plasgenomicsutils variant_spacing --input calls.bcf\n"
               "    plasgenomicsutils variant_spacing --input calls.bcf \\\n"
               "        --locus Pf3D7_07_v3:886682-891682,Pf3D7_07_v3:889636-889836\n"
               "    plasgenomicsutils variant_spacing --input calls.bcf --bed panel.bed \\\n"
               "        --out spacing.tsv\n\n"
               "Gaps are measured within a chromosome; the step from the end of one to the\n"
               "start of the next is never counted. The `all` row pools every chromosome's\n"
               "gaps -- it is not the mean of the rows above it.\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True, help="VCF/BCF (indexed, if regions are used)")
    p.add_argument("--out", default=None, help="Write here instead of stdout")

    r = p.add_argument_group("restricting to part of the genome")
    r.add_argument("--locus", default=None,
                   help="Comma-separated regions, e.g. "
                        "'Pf3D7_07_v3:886682-891682,Pf3D7_07_v3:889636-889836'. These are "
                        "**1-based inclusive**, as bcftools -r reads them. Overlapping "
                        "regions are unioned for the density denominator, so a nested one "
                        "does not deflate it.")
    r.add_argument("--bed", default=None,
                   help="BED of regions instead of --locus. **0-based half-open**, as any "
                        "BED is -- the opposite convention to --locus, which is bcftools' "
                        "own split and kept here so each reads the way its format does.")

    d = p.add_argument_group("density")
    d.add_argument("--reference", default=DEFAULT_REFERENCE,
                   help=f"Reference for the genetic-map rate (default: {DEFAULT_REFERENCE}; "
                        f"available: {', '.join(available_references())})")
    d.add_argument("--bp-per-cm", type=float, default=None,
                   help="Constant genetic-map rate (bp/cM), overriding the reference's "
                        "(Pf3D7 = 15000)")
    return p


def parse_args_variant_spacing():
    return get_parser_variant_spacing().parse_args()


def variant_spacing():
    args = parse_args_variant_spacing()
    bp_per_cm = (args.bp_per_cm if args.bp_per_cm is not None
                 else get_reference(args.reference).bp_per_cm)
    rows = _spacing(args.input, locus=args.locus, bed=args.bed, bp_per_cm=bp_per_cm)
    if not rows:
        raise SystemExit("ERROR: this callset has no variants at all")
    write_table(rows, args.out)
    if args.out:
        print(f"  -> {args.out}")


if __name__ == "__main__":
    variant_spacing()
