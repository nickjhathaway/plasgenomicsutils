#!/usr/bin/env python
"""Compute global and per-region alternate allele frequencies from a BCF/VCF."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ...lib.ibd_freqs import compute_allele_freqs as _compute
from ...utils.small_utils import Utils


def get_parser_compute_allele_freqs() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils compute_allele_freqs",
        description="Compute global and per-region allele frequencies from a BCF/VCF",
    )
    p.add_argument("--bcf", required=True, help="Input BCF/VCF")
    p.add_argument("--meta", required=True,
                   help="Table with at least columns: sample, <region-col>")
    p.add_argument("--region-col", default="region",
                   help="Column in --meta to use as region (default: region)")
    p.add_argument("--zero-based", action="store_true",
                   help="Emit 0-based snp_ids (chr:pos-1) to match build_ibd_matrix labels")
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

    coord = "0-based" if args.zero_based else "1-based (VCF native)"
    print(f"SNP ID coordinate system: {coord}")

    sample_to_region = (
        meta.dropna(subset=[args.region_col])
            .set_index("sample")[args.region_col]
            .astype(str)
            .to_dict()
    )

    print("\nComputing global + per-region AF (single pass)...")
    global_df, region_df = _compute(
        args.bcf, sample_to_region=sample_to_region, zero_based=args.zero_based
    )

    global_out = outdir / "allele_freqs.tsv.gz"
    Utils.write_tsv_gz(global_df, str(global_out))
    print(f"  -> {global_out}  ({len(global_df):,} SNPs)")
    print(f"  Example SNP IDs: {global_df['snp_id'].head(3).tolist()}")

    region_out = outdir / "region_allele_freqs.tsv.gz"
    Utils.write_tsv_gz(region_df, str(region_out))
    n_regions = region_df["region"].nunique() if len(region_df) else 0
    print(f"  -> {region_out}  ({len(region_df):,} rows, {n_regions} regions)")


if __name__ == "__main__":
    compute_allele_freqs()
