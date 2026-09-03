#!/usr/bin/env python
"""Compute alternate allele frequencies from a BCF/VCF, overall and optionally per group."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ...lib.ibd_freqs import compute_allele_freqs as _compute
from ...lib.intervals import SNP_COORD_SYSTEM
from ...utils.small_utils import Utils


def get_parser_compute_allele_freqs() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils compute_allele_freqs",
        description="Compute allele frequencies from a BCF/VCF: one table over every "
                    "sample, plus a per-group table when --meta is given.",
    )
    p.add_argument("--bcf", required=True, help="Input BCF/VCF")
    p.add_argument("--meta", default=None,
                   help="Table with at least columns: sample, <group-col>. Optional -- "
                        "without it only the whole-file table is written, which is the "
                        "quickest way to look at what frequencies a callset actually holds.")
    p.add_argument("--group-col", default="group",
                   help="Column in --meta to use as group (default: group); needs --meta")
    p.add_argument("--per-alt", action="store_true",
                   help="One row per (SNP, ALT) instead of one per SNP, with alt and "
                        "alt_index columns. At a biallelic site the two agree; at a "
                        "multiallelic one the default collapses every ALT together and "
                        "this separates them. Leave it off for tables the selection "
                        "statistic joins on, which expect one row per SNP.")
    p.add_argument("--ad-min-reads", type=int, default=2,
                   help="Reads supporting an allele for prevalence_ad to count the sample "
                        "as carrying it (default: 2, as filter_ad_regenotype).")
    p.add_argument("--ad-min-freq", type=float, default=0.01,
                   help="Within-sample fraction for the same (default: 0.01). Raise it to "
                        "ignore the lowest-frequency clones.")
    p.add_argument("--no-weighted", action="store_true",
                   help="Skip af_weighted (the mean within-sample frequency, from "
                        "FORMAT/AD). Saves an array read per record; the column is NaN "
                        "anyway where AD is absent.")
    p.add_argument("--with-pos-vcf", action="store_true",
                   help="Also emit a pos_vcf column (1-based VCF position) for looking "
                        "variants up by eye; off by default to keep the table small")
    p.add_argument("--output", default=".", help="Output directory (default: cwd)")
    p.add_argument("--meta-file-separator", default="tab",
                   help="Separator of --meta: 'tab', 'comma', or a literal char (default: tab)")
    return p


def parse_args_compute_allele_freqs():
    return get_parser_compute_allele_freqs().parse_args()


def compute_allele_freqs():
    args = parse_args_compute_allele_freqs()

    outdir = Utils.ensure_dir(args.output)
    print(f"SNP ID coordinate system: {SNP_COORD_SYSTEM} (chr:pos0)")

    sample_to_group = None
    if args.meta:
        sep = Utils.resolve_delim(args.meta_file_separator)
        # through read_meta, not pd.read_csv: it is what makes the sample ids strings, and
        # an all-numeric cohort read as int64 matches none of the VCF's sample names
        meta = Utils.read_meta(args.meta, sep=sep, wants=("sample", args.group_col))
        for col in ("sample", args.group_col):
            Utils.resolve_column(meta.columns, col, source=f"metadata ({args.meta})")
        sample_to_group = (
            meta.dropna(subset=[args.group_col])
                .set_index("sample")[args.group_col]
                .astype(str)
                .to_dict()
        )
    elif args.group_col != "group":
        raise SystemExit("ERROR: --group-col needs --meta")

    print("\nComputing AF over every sample"
          + (" + per group" if sample_to_group else "") + " (single pass)...")
    global_df, group_df = _compute(
        args.bcf, sample_to_group=sample_to_group, with_pos_vcf=args.with_pos_vcf,
        per_alt=args.per_alt, weighted=not args.no_weighted,
        ad_min_reads=args.ad_min_reads, ad_min_freq=args.ad_min_freq,
    )

    global_out = outdir / "allele_freqs.tsv.gz"
    Utils.write_tsv_gz(global_df, str(global_out),
                       header_comment=f"snp_coord_system={SNP_COORD_SYSTEM}")
    unit = "rows (SNP x ALT)" if args.per_alt else "SNPs"
    print(f"  -> {global_out}  ({len(global_df):,} {unit})")
    print(f"  Example SNP IDs: {global_df['snp_id'].head(3).tolist()}")

    if sample_to_group is None:
        # no grouping asked for, so no empty second file to explain later
        print("  no --meta, so no per-group table written")
        return

    group_out = outdir / "group_allele_freqs.tsv.gz"
    Utils.write_tsv_gz(group_df, str(group_out),
                       header_comment=f"snp_coord_system={SNP_COORD_SYSTEM}")
    n_groups = group_df["group"].nunique() if len(group_df) else 0
    print(f"  -> {group_out}  ({len(group_df):,} rows, {n_groups} groups)")


if __name__ == "__main__":
    compute_allele_freqs()
