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

    print("=== harmonize_bcf ===", file=sys.stderr)
    print("  Input files:", file=sys.stderr)
    for fpath, out_path in zip(args.files, out_paths):
        print(f"    {fpath}  ->  {out_path}", file=sys.stderr)
    print(f"  min_ad={args.min_ad} (non-inclusive)  min_af={args.min_af} (non-inclusive)  "
          f"het_min_af={args.het_min_af}", file=sys.stderr)

    print("\n=== Step 1: cleaning spurious low-level ALTs, building the ALT union ===",
          file=sys.stderr)
    if drop_indels:
        print("  dropping indel-context records from all inputs (--keep-indels to disable)",
              file=sys.stderr)
    union, dup_positions, ambiguous, p1 = H.accumulate_union(
        args.files, args.min_ad, args.min_af, args.het_min_af, drop_indels=drop_indels)
    for fpath in args.files:
        st = p1["per_file"][fpath]
        print(f"\n  [{os.path.basename(fpath)}]", file=sys.stderr)
        print(f"    Records read:                     {st['processed']:,}", file=sys.stderr)
        if drop_indels and st["indel_context"]:
            print(f"    Indel-context records dropped:    {st['indel_context']:,}",
                  file=sys.stderr)
        print(f"    Sites kept:                       {st['sites']:,}", file=sys.stderr)
        print(f"    ALT alleles zeroed and removed:   {st['alts_removed']:,}", file=sys.stderr)
        print(f"    Records reduced to ref-only:      {st['reduced_to_ref_only']:,} "
              f"(kept for the union)", file=sys.stderr)

    print(f"\n  Total union sites:                  {p1['union_sites']:,}", file=sys.stderr)
    print(f"  Sites with at least one ALT:        {p1['union_with_alts']:,}", file=sys.stderr)
    print(f"  Sites dropped (ref-only everywhere):{p1['union_dropped']:,}", file=sys.stderr)
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

    print("\n=== Step 2: harmonizing each file to the union ===", file=sys.stderr)
    absent = {}
    for fpath, out_path in zip(args.files, out_paths):
        print(f"\n  [{os.path.basename(fpath)}] -> [{os.path.basename(out_path)}]",
              file=sys.stderr)
        if not need_bcftools:
            st = H.harmonize_file(fpath, out_path, union, args.min_ad, args.min_af,
                                  args.het_min_af, drop_indels=drop_indels)
        else:
            tmp_vcf = out_path + ".tmp.vcf"
            st = H.harmonize_file(fpath, tmp_vcf, union, args.min_ad, args.min_af,
                                  args.het_min_af, drop_indels=drop_indels)
            # annotate strips stale FORMAT and writes the requested output format
            sh(f"bcftools annotate {strip_arg} -O{args.output_format} -o {q(out_path)} {q(tmp_vcf)}",
               tools=("bcftools",))
            os.remove(tmp_vcf)
        absent[fpath] = st["absent"]
        print(f"    Records written:                  {st['written']:,}", file=sys.stderr)
        print(f"    Records harmonized (ALTs added):  {st['alts_added']:,}", file=sys.stderr)
        print(f"    Records dropped (ref-only):       {st['dropped_ref_only']:,}",
              file=sys.stderr)
        print(f"    Union sites absent from the file: {st['absent']:,}", file=sys.stderr)

    # The number above is the one that bites later, so spell out what it means rather than
    # leaving it as a statistic: harmonizing settles the ALT sets, not which sites each file
    # holds, so bcftools merge fills the gaps with missing genotypes AND missing FORMAT/AD.
    # hmmibd-rs reads AD as an integer and stops at the first missing one
    # ("NumericaValueEmptyInt"), so this is worth acting on before the merge, not after.
    if any(absent.values()):
        worst = max(absent.values())
        print(f"\n  NOTE: up to {worst:,} union site(s) are absent from an input, so those "
              f"samples get\n        missing GT *and* missing FORMAT/AD after `bcftools merge`. "
              f"Tools that read AD as\n        an integer (hmmibd-rs) fail on that. To keep only "
              f"sites every cohort called:\n"
              f"          bcftools view -e 'FMT/AD=\".\"' merged.bcf -Ob -o merged.nomiss.bcf",
              file=sys.stderr)

    print("\n=== Done ===", file=sys.stderr)


if __name__ == "__main__":
    harmonize_bcf()
