#!/usr/bin/env python
"""Linkage-disequilibrium decay: mean r-squared against SNP-pair distance, per group."""

from __future__ import annotations

import argparse

import pandas as pd

from ...lib.ld import (
    LD_MAX_DIST,
    LD_MAX_SNPS,
    LD_MIN_MAF,
    ld_decay as compute_ld_decay,
    read_dosages,
)
from ...utils.small_utils import Utils


def get_parser_ld_decay() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils ld_decay",
        description="Mean r-squared between SNP pairs, binned by the distance between "
                    "them, for each metadata group. How fast it falls says how freely the "
                    "population recombines",
    )
    p.add_argument("--vcf", required=True, help="VCF/BCF (the full, unpruned callset -- "
                                                "LD-pruning removes what this measures)")
    p.add_argument("--meta", default=None,
                   help="Sample metadata; without it every sample is pooled")
    p.add_argument("--group-col", default="region",
                   help="Metadata column defining the groups (default: %(default)s)")
    p.add_argument("--region", action="append", default=None, metavar="REG",
                   help="chrom or chrom:start-end; repeatable (needs an index)")
    p.add_argument("--max-dist", type=int, default=LD_MAX_DIST,
                   help="Longest SNP-pair separation to count, bp (default: %(default)s)")
    p.add_argument("--bins", type=int, default=20,
                   help="Number of distance bins (default: %(default)s)")
    p.add_argument("--maf", type=float, default=LD_MIN_MAF,
                   help="Minor-allele-frequency floor within each group; rare alleles "
                        "give noisy, systematically low r-squared (default: %(default)s)")
    p.add_argument("--max-snps", type=int, default=LD_MAX_SNPS,
                   help="SNPs kept per chromosome before the pairwise scan, evenly "
                        "spaced. The scan is quadratic in SNPs inside one window "
                        "(default: %(default)s)")
    p.add_argument("--min-samples", type=int, default=4,
                   help="Skip groups smaller than this (default: %(default)s)")
    p.add_argument("--min-depth", type=int, default=0,
                   help="Treat a genotype backed by fewer reads as uncalled "
                        "(default: %(default)s, 0 disables)")
    p.add_argument("--het", choices=("missing", "dosage"), default="missing",
                   help="How to read a heterozygous call. The parasite is haploid, so "
                        "'missing' treats it as a mixed infection (default: %(default)s)")
    p.add_argument("--output", required=True, help="Binned r-squared TSV(.gz)")
    p.add_argument("--half-decay-output", default=None,
                   help="Per-group half-decay distance TSV")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    return p


def parse_args_ld_decay():
    return get_parser_ld_decay().parse_args()


def ld_decay():
    args = parse_args_ld_decay()
    Utils.output_file_check(args.output, args.overwrite)
    if args.half_decay_output:
        Utils.output_file_check(args.half_decay_output, args.overwrite)

    samples = groups = None
    if args.meta:
        meta = Utils.read_table(args.meta)
        if args.group_col not in meta.columns:
            raise SystemExit(f"--meta has no '{args.group_col}' column")
        meta = meta.dropna(subset=["sample", args.group_col])
        samples = meta["sample"].astype(str).tolist()

    gn, chrom, pos, names = read_dosages(args.vcf, samples=samples, regions=args.region,
                                         het=args.het, min_depth=args.min_depth)
    print(f"[info] {gn.shape[0]:,} variants x {gn.shape[1]} samples")
    if args.meta:
        groups = meta.set_index("sample")[args.group_col].reindex(names).to_numpy()
        missing = int(pd.isna(groups).sum())
        if missing:
            print(f"[info] {missing} sample(s) have no {args.group_col} and are pooled out")

    df, half = compute_ld_decay(gn, chrom, pos, groups=groups, max_dist=args.max_dist,
                                bins=args.bins, maf=args.maf, max_snps=args.max_snps,
                                min_samples=args.min_samples)
    if df.empty:
        raise SystemExit("no group produced any SNP pairs; loosen --maf or --min-samples")

    header = (f"#max_dist={args.max_dist}\t#max_snps={args.max_snps}\t#maf={args.maf}")
    Utils.write_tsv_gz(df, args.output, header_comment=header)
    print(f"[done] wrote {args.output} ({len(df)} rows)")
    for r in half.itertuples(index=False):
        d = "still above half at max_dist" if pd.isna(r.half_decay_bp) \
            else f"{r.half_decay_bp:,.0f} bp"
        print(f"  {r.group}: half-decay {d}")
    if args.half_decay_output:
        half.to_csv(args.half_decay_output, sep="\t", index=False)
        print(f"[done] wrote {args.half_decay_output}")


if __name__ == "__main__":
    ld_decay()
