#!/usr/bin/env python
"""Regions almost no sample covers -- candidate sWGA amplification dropouts."""

from __future__ import annotations

import argparse

import pandas as pd

from ...lib.assets import resolve_bed
from ...lib.coverage import (
    DROPOUT_MIN_DEPTH,
    DROPOUT_MIN_FRAC_SAMPLES,
    annotate_regions,
    dropout_regions as find_dropouts,
    load_bed,
)
from ...utils.small_utils import Utils


def get_parser_dropout_regions() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils coverage_dropout_regions",
        description="Find windows that are below depth in nearly every sample and merge "
                    "them into regions. Selective whole-genome amplification does not "
                    "amplify uniformly, and a region no sample amplifies reads as "
                    "invariant rather than as missing -- so flag it instead",
    )
    p.add_argument("--windows", required=True, action="append", metavar="TSV",
                   help="Per-window depth table from coverage_depth_stats "
                        "(--windows-output); repeatable to combine batches")
    p.add_argument("--min-depth", type=float, default=DROPOUT_MIN_DEPTH,
                   help="A sample counts as uncovered in a window below this mean depth "
                        "(default: %(default)s)")
    p.add_argument("--min-frac-samples", type=float, default=DROPOUT_MIN_FRAC_SAMPLES,
                   help="Flag a window when at least this fraction of samples are "
                        "uncovered (default: %(default)s)")
    p.add_argument("--merge-gap", type=int, default=0,
                   help="Join flagged regions separated by at most this many bp "
                        "(default: %(default)s)")
    p.add_argument("--min-length", type=int, default=0,
                   help="Drop merged regions shorter than this many bp "
                        "(default: %(default)s)")
    p.add_argument("--regions", default=None, metavar="BED",
                   help="Keep only windows inside this BED (accepts "
                        "builtin:pf3d7_core_regions). Subtelomeric dropout is expected, "
                        "so restricting to the core is what makes the output actionable")
    p.add_argument("--genes", default=None,
                   help="Gene intervals (name, chr, start, end) to annotate each region "
                        "with the genes it overlaps")
    p.add_argument("--output", required=True, help="Merged dropout regions TSV(.gz)")
    p.add_argument("--bed-output", default=None,
                   help="Also write the regions as a plain BED, ready to mask with")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    return p


def parse_args_dropout_regions():
    return get_parser_dropout_regions().parse_args()


def _in_regions(windows, regions):
    keep = []
    for chrom, sub in windows.groupby("chrom", sort=False):
        iv = regions.get(chrom, [])
        if not iv:
            continue
        mask = pd.Series(False, index=sub.index)
        for start, end in iv:
            mask |= (sub["start"] < end) & (sub["end"] > start)
        keep.append(sub[mask])
    return pd.concat(keep, ignore_index=True) if keep else windows.iloc[0:0]


def dropout_regions():
    args = parse_args_dropout_regions()
    Utils.output_file_check(args.output, args.overwrite)
    if args.bed_output:
        Utils.output_file_check(args.bed_output, args.overwrite)

    windows = pd.concat([Utils.read_table(w) for w in args.windows], ignore_index=True)
    need = {"sample", "chrom", "start", "end", "mean_depth"}
    missing = need - set(windows.columns)
    if missing:
        raise SystemExit(f"--windows is missing column(s): {', '.join(sorted(missing))}")
    n_samples = windows["sample"].nunique()
    if args.regions:
        windows = _in_regions(windows, load_bed(resolve_bed(args.regions)))

    regions = find_dropouts(windows, min_depth=args.min_depth,
                            min_frac_samples=args.min_frac_samples,
                            merge_gap=args.merge_gap, min_length=args.min_length)
    if args.genes and not regions.empty:
        genes = Utils.read_table(args.genes)
        if "chr" not in genes.columns and "chrom" in genes.columns:
            genes = genes.rename(columns={"chrom": "chr"})
        regions = annotate_regions(regions, genes)

    Utils.write_tsv_gz(regions, args.output)
    total_bp = int(regions["length"].sum()) if not regions.empty else 0
    print(f"[done] wrote {args.output}: {len(regions)} region(s), {total_bp:,} bp "
          f"below {args.min_depth}x in >= {args.min_frac_samples:.0%} of {n_samples} samples")
    if args.bed_output:
        cols = ["chrom", "start", "end"]
        regions[cols].to_csv(args.bed_output, sep="\t", header=False, index=False)
        print(f"[done] wrote {args.bed_output}")


if __name__ == "__main__":
    dropout_regions()
