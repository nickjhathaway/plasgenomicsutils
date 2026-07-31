#!/usr/bin/env python
"""Harmonize the ALT sets of separately-called cohorts so bcftools can merge them."""

from __future__ import annotations

import argparse
import os
import sys

from ...lib import harmonize as H
from ...lib.bcftools import q, require, sh


def _output_path(stub: str, input_path: str, fmt: str) -> str:
    basename = os.path.basename(input_path)
    for ext in (".bcf", ".vcf.gz", ".vcf"):
        if basename.endswith(ext):
            basename = basename[: -len(ext)]
            break
    return f"{stub}_{basename}{H.OUTPUT_EXT[fmt]}"


def get_parser_harmonize_bcf() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils harmonize_bcf",
        description="Clean spurious ALTs and harmonize allele sets across files for merging.",
    )
    p.add_argument("--files", required=True, nargs="+",
                   help="Two or more coordinate-sorted BCF/VCF input files")
    p.add_argument("--stub", required=True,
                   help="Output stub; each output is STUB_<input_basename>.vcf")
    p.add_argument("--min-ad", type=int, default=3,
                   help="Zero ALT AD below this read count (default: 3)")
    p.add_argument("--min-af", type=float, default=0.005,
                   help="Zero ALT AD below this within-sample frequency (default: 0.005)")
    p.add_argument("--het-min-af", type=float, default=0.2,
                   help="Minimum minor-allele frequency to call a heterozygote (default: 0.2)")
    p.add_argument("--output-format", choices=["v", "z", "b"], default="v",
                   help="Output format: v=VCF (default), z=bgzipped VCF, b=BCF. "
                        "z/b are produced by writing VCF then converting with bcftools; "
                        "harmonize never writes BCF directly, because pysam leaves a "
                        "binary AD/allele mismatch when a record's alleles are reduced.")
    p.add_argument("--keep-indels", action="store_true",
                   help="Keep indel-context records. By default they are dropped from every "
                        "input — including no-ALT 'REF>.' records carrying the INDEL flag, "
                        "which `bcftools view --exclude-type indels` does not remove.")
    return p


def parse_args_harmonize_bcf():
    return get_parser_harmonize_bcf().parse_args()


def harmonize_bcf():
    args = parse_args_harmonize_bcf()
    if len(args.files) < 2:
        raise SystemExit("ERROR: at least 2 input files are required")

    out_paths = [_output_path(args.stub, f, args.output_format) for f in args.files]
    drop_indels = not args.keep_indels
    # FORMAT fields (e.g. PL) that go stale once alleles are reshaped; harmonize
    # maintains AD but cannot recompute these, so drop them or bcftools merge fails.
    stale_fmt = sorted(set().union(*(set(H.stale_format_fields(f)) for f in args.files)))
    need_bcftools = args.output_format != "v" or bool(stale_fmt)
    if need_bcftools:
        require("bcftools")  # for the strip and/or VCF -> z/b conversion
    if stale_fmt:
        print(f"  stripping allele-dependent FORMAT field(s) invalidated by harmonizing: "
              f"{', '.join(stale_fmt)}", file=sys.stderr)

    print("=== Pass 1: cleaning + building ALT union ===", file=sys.stderr)
    if drop_indels:
        print("  dropping indel-context records from all inputs (--keep-indels to disable)",
              file=sys.stderr)
    union, dup_positions, ambiguous = H.accumulate_union(
        args.files, args.min_ad, args.min_af, args.het_min_af, drop_indels=drop_indels)
    n_with_alts = sum(1 for v in union.values() if len(v) > 1)
    print(f"  union sites with >= 1 ALT: {n_with_alts:,}", file=sys.stderr)
    if dup_positions:
        print(f"  NOTE: {len(dup_positions)} duplicate (chrom,pos) record(s) collapsed; kept the "
              f"record with the most real ALT alleles (drops overlapping no-ALT/indel records).",
              file=sys.stderr)
        for fpath, chrom, pos in sorted(dup_positions)[:10]:
            print(f"    {chrom}:{pos}  in {os.path.basename(fpath)}", file=sys.stderr)
    if ambiguous:
        print(f"  WARNING: {len(ambiguous)} position(s) had >1 record carrying real ALTs "
              f"(un-normalized). Run `bcftools norm -m -any` on inputs first for a correct union.",
              file=sys.stderr)
        for fpath, chrom, pos in sorted(ambiguous)[:10]:
            print(f"    {chrom}:{pos}  in {os.path.basename(fpath)}", file=sys.stderr)

    strip_arg = ("-x " + q("FORMAT/" + ",FORMAT/".join(stale_fmt))) if stale_fmt else ""

    print("=== Pass 2: harmonizing each file to the union ===", file=sys.stderr)
    for fpath, out_path in zip(args.files, out_paths):
        print(f"  {fpath} -> {out_path}", file=sys.stderr)
        if not need_bcftools:
            H.harmonize_file(fpath, out_path, union, args.min_ad, args.min_af,
                             args.het_min_af, drop_indels=drop_indels)
        else:
            tmp_vcf = out_path + ".tmp.vcf"
            H.harmonize_file(fpath, tmp_vcf, union, args.min_ad, args.min_af,
                             args.het_min_af, drop_indels=drop_indels)
            # annotate strips stale FORMAT and writes the requested output format
            sh(f"bcftools annotate {strip_arg} -O{args.output_format} -o {q(out_path)} {q(tmp_vcf)}",
               tools=("bcftools",))
            os.remove(tmp_vcf)

    print("=== Done ===", file=sys.stderr)


if __name__ == "__main__":
    harmonize_bcf()
