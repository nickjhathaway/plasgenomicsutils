#!/usr/bin/env python
"""Per-sample singleton counts, with an outlier flag for sample QC."""

from __future__ import annotations

import argparse

from ...lib.singletons import (DEFAULT_DUPLICATE_FRAC, DEFAULT_MAD_CUTOFF,
                               count_singletons, flag_outliers)
from ...utils.small_utils import Utils


def get_parser_singleton_counts() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils singleton_counts",
        description="Count, per sample, the variants where it is the only non-reference "
                    "carrier. A sample carrying far more private variants than the rest "
                    "of the cohort is usually contaminated, mixed-species or "
                    "mis-aligned rather than interesting",
    )
    p.add_argument("--vcf", required=True, help="VCF/BCF (indexed if --region is used)")
    p.add_argument("--samples", default=None,
                   help="File of sample names, one per line, or a comma-separated list; "
                        "singleton status is judged within the set analysed")
    p.add_argument("--region", action="append", default=None, metavar="REG",
                   help="chrom or chrom:start-end; repeatable")
    p.add_argument("--max-missing-frac", type=float, default=0.2,
                   help="Skip variants with more missing genotypes than this; a variant "
                        "called in few samples makes a private allele trivially likely "
                        "(default: %(default)s)")
    p.add_argument("--min-depth", type=int, default=5,
                   help="Treat a genotype backed by fewer than this many reads as "
                        "uncalled. Callsets built from gVCFs emit 0/0 at sites with no "
                        "reads at all rather than ./., which would otherwise count as a "
                        "confident reference call; depth is read from AD (summed) or DP "
                        "(default: %(default)s, 0 disables)")
    p.add_argument("--mad-cutoff", type=float, default=DEFAULT_MAD_CUTOFF,
                   help="Flag samples this many MADs above the cohort median rate "
                        "(default: %(default)s)")
    p.add_argument("--duplicate-frac", type=float, default=DEFAULT_DUPLICATE_FRAC,
                   help="Call two samples near-identical when this share of one's "
                        "doubletons (variants carried by exactly two samples) are shared "
                        "with a single partner. Such a sample has almost no singletons, "
                        "because its private variants became the pair's. Strict on "
                        "purpose: merely related samples reach ~0.6-0.8 "
                        "(default: %(default)s)")
    p.add_argument("--output", required=True, help="Per-sample TSV(.gz)")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing output")
    return p


def parse_args_singleton_counts():
    return get_parser_singleton_counts().parse_args()


def _sample_list(value):
    if not value:
        return None
    from pathlib import Path
    if Path(value).exists():
        return [s.strip() for s in Path(value).read_text().splitlines() if s.strip()]
    return [s.strip() for s in value.split(",") if s.strip()]


def singleton_counts():
    args = parse_args_singleton_counts()
    Utils.output_file_check(args.output, args.overwrite)

    df, n_variants = count_singletons(args.vcf, samples=_sample_list(args.samples),
                                      regions=args.region,
                                      max_missing_frac=args.max_missing_frac,
                                      min_depth=args.min_depth)
    n_low = df.attrs.get("n_low_depth", 0)
    df = flag_outliers(df, mad_cutoff=args.mad_cutoff,
                       duplicate_frac=args.duplicate_frac)
    df = df.sort_values("singleton_rate", ascending=False).reset_index(drop=True)
    Utils.write_tsv_gz(df, args.output)

    n_out = int(df["outlier"].sum())
    dups = df[df["flag"].str.contains("near-identical", na=False)]
    print(f"[done] wrote {args.output}: {len(df)} samples over {n_variants:,} variants; "
          f"{int(df['n_singleton'].sum()):,} singletons, {n_out} outlier(s)")
    if n_low:
        print(f"       {n_low:,} genotype(s) read as uncalled below {args.min_depth} reads")
    for r in df[df["flag"] != ""].itertuples(index=False):
        print(f"  {r.sample}\t{r.n_singleton:,} singletons\t"
              f"{r.singleton_rate:.2f}/1000 called\t{r.mad_score:+.1f} MADs\t{r.flag}")
    if len(dups):
        print(f"       {len(dups)} sample(s) are near-identical to another -- the same "
              f"parasite sequenced twice, or one clone in two hosts;")
        print(f"       check IBD and the collection records before dropping either")


if __name__ == "__main__":
    singleton_counts()
