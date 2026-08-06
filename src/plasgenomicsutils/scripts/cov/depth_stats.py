#!/usr/bin/env python
"""Per-sample sequencing depth summaries from BAM/CRAM."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ...lib.assets import resolve_bed
from ...lib.coverage import (
    DEFAULT_THRESHOLDS,
    DEFAULT_WINDOW,
    load_bed,
    mosdepth_available,
    resolve_engine,
    sample_coverage,
)
from ...utils.small_utils import Utils


def get_parser_depth_stats() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils coverage_depth_stats",
        description="Depth of coverage per sample: mean/median/SD/quartiles and breadth "
                    "at a set of thresholds, per chromosome and genome-wide, optionally "
                    "restricted to a BED (e.g. the core genome)",
    )
    p.add_argument("--bam", action="append", default=None, metavar="PATH",
                   help="Indexed BAM/CRAM; repeatable")
    p.add_argument("--bam-list", default=None,
                   help="File of BAM/CRAM paths, one per line; optional second "
                        "tab-separated column gives the sample name")
    p.add_argument("--regions", default=None, metavar="BED",
                   help="Count only bases inside this BED (accepts "
                        "builtin:pf3d7_core_regions). Subtelomeric and hypervariable "
                        "regions otherwise drag every statistic down")
    p.add_argument("--chrom", action="append", default=None, metavar="NAME",
                   help="Restrict to this contig; repeatable")
    p.add_argument("--thresholds", default=",".join(str(t) for t in DEFAULT_THRESHOLDS),
                   help="Comma-separated depth cutoffs to report breadth for "
                        "(default: %(default)s)")
    p.add_argument("--window", type=int, default=DEFAULT_WINDOW,
                   help="Window size in bp for the per-window means that feed "
                        "coverage_dropout_regions (default: %(default)s)")
    p.add_argument("--engine", choices=("auto", "pysam", "mosdepth"), default="auto",
                   help="Depth engine; auto uses mosdepth when it is on PATH and falls "
                        "back to pysam (default: %(default)s)")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="Ignore reads below this mapping quality (default: %(default)s)")
    p.add_argument("--min-baseq", type=int, default=0,
                   help="Ignore bases below this base quality; pysam engine only "
                        "(default: %(default)s)")
    p.add_argument("--threads", type=int, default=1,
                   help="Threads for the mosdepth engine (default: %(default)s)")
    p.add_argument("--jobs", type=int, default=1,
                   help="Samples to process in parallel. The pysam engine reads roughly "
                        "90 kb of genome per second per core, so a whole-genome cohort is "
                        "worth spreading out (default: %(default)s)")
    p.add_argument("--reference", default=None,
                   help="Reference FASTA, required for CRAM")
    p.add_argument("--output", required=True,
                   help="Per-sample, per-chromosome summary TSV(.gz)")
    p.add_argument("--windows-output", default=None,
                   help="Per-window mean depth TSV(.gz); the input to "
                        "coverage_dropout_regions")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    return p


def parse_args_depth_stats():
    return get_parser_depth_stats().parse_args()


def _bam_entries(args):
    entries: list[tuple[str, str | None]] = []
    for b in args.bam or []:
        entries.append((b, None))
    if args.bam_list:
        with open(args.bam_list) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                entries.append((parts[0], parts[1] if len(parts) > 1 else None))
    if not entries:
        raise SystemExit("give at least one --bam or a --bam-list")
    missing = [b for b, _ in entries if not Path(b).exists()]
    if missing:
        raise SystemExit(f"BAM/CRAM not found: {', '.join(missing[:5])}")
    return entries


def depth_stats():
    args = parse_args_depth_stats()
    Utils.output_file_check(args.output, args.overwrite)
    if args.windows_output:
        Utils.output_file_check(args.windows_output, args.overwrite)

    entries = _bam_entries(args)
    thresholds = tuple(int(t) for t in args.thresholds.split(",") if t.strip())
    regions = load_bed(resolve_bed(args.regions)) if args.regions else None
    engine = resolve_engine(args.engine, args.min_baseq)
    if args.engine == "auto" and engine == "pysam" and not mosdepth_available():
        print("[info] mosdepth not found; using the pysam engine", flush=True)

    kwargs = dict(regions=regions, window=args.window, thresholds=thresholds,
                  chroms=args.chrom, engine=engine, min_mapq=args.min_mapq,
                  min_baseq=args.min_baseq, threads=args.threads,
                  reference=args.reference)

    def report(per_chrom):
        g = next(r for r in per_chrom if r["chrom"] == "genome")
        print(f"[info] {g['sample']}: mean {g['mean']:.1f}x, median {g['median']:.0f}x, "
              f"{g[f'pct_ge_{thresholds[0]}x']:.1f}% >= {thresholds[0]}x", flush=True)

    rows, windows = [], []
    if args.jobs > 1 and len(entries) > 1:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(sample_coverage, path, sample=name, **kwargs)
                       for path, name in entries]
            for fut in futures:
                per_chrom, win = fut.result()
                rows.extend(per_chrom)
                if args.windows_output:
                    windows.extend(win)
                report(per_chrom)
    else:
        for path, name in entries:
            per_chrom, win = sample_coverage(path, sample=name, **kwargs)
            rows.extend(per_chrom)
            if args.windows_output:
                windows.extend(win)
            report(per_chrom)

    Utils.write_tsv_gz(pd.DataFrame(rows), args.output)
    print(f"[done] wrote {args.output} ({len(rows)} rows)")
    if args.windows_output:
        Utils.write_tsv_gz(pd.DataFrame(windows), args.windows_output)
        print(f"[done] wrote {args.windows_output} ({len(windows)} rows)")


if __name__ == "__main__":
    depth_stats()
