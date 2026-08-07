#!/usr/bin/env python
"""Compute global and per-group alternate allele frequencies from a BCF/VCF."""

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
        description="Compute global and per-group allele frequencies from a BCF/VCF",
    )
    p.add_argument("--bcf", required=True, help="Input BCF/VCF")
    p.add_argument("--meta", required=True,
                   help="Table with at least columns: sample, <group-col>")
    p.add_argument("--group-col", default="group",
                   help="Column in --meta to use as group (default: group)")
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

    sep = Utils.resolve_delim(args.meta_file_separator)
    meta = pd.read_csv(args.meta, sep=sep)
    outdir = Utils.ensure_dir(args.output)

    print(f"SNP ID coordinate system: {SNP_COORD_SYSTEM} (chr:pos0)")

    sample_to_group = (
        meta.dropna(subset=[args.group_col])
            .set_index("sample")[args.group_col]
            .astype(str)
            .to_dict()
    )

    print("\nComputing global + per-group AF (single pass)...")
    global_df, group_df = _compute(
        args.bcf, sample_to_group=sample_to_group, with_pos_vcf=args.with_pos_vcf
    )

    global_out = outdir / "allele_freqs.tsv.gz"
    Utils.write_tsv_gz(global_df, str(global_out),
                       header_comment=f"snp_coord_system={SNP_COORD_SYSTEM}")
    print(f"  -> {global_out}  ({len(global_df):,} SNPs)")
    print(f"  Example SNP IDs: {global_df['snp_id'].head(3).tolist()}")

    group_out = outdir / "group_allele_freqs.tsv.gz"
    Utils.write_tsv_gz(group_df, str(group_out),
                       header_comment=f"snp_coord_system={SNP_COORD_SYSTEM}")
    n_groups = group_df["group"].nunique() if len(group_df) else 0
    print(f"  -> {group_out}  ({len(group_df):,} rows, {n_groups} groups)")


if __name__ == "__main__":
    compute_allele_freqs()
