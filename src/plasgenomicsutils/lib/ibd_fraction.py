"""Per-pair IBD fraction and SNP density from hmmibd-rs blocks.

Chromosome lengths and the genetic-map rate come from
:mod:`plasgenomicsutils.lib.reference` (selected by ``--reference``), so other
species can be added without touching this math.

Callable-genome denominator: Pf sub-telomeres cannot be reliably SNP-called, so
IBD can never be detected there. The denominator is the callable span
(per-chromosome last-SNP minus first-SNP, summed), not full chromosome length.
Under a constant map rate the rate cancels: f = total_IBD_bp / callable_bp.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .ibd_matrix import IBD_MIN_BLOCK_KB, IBD_MIN_BLOCK_SNP, filter_ibd_blocks, read_blocks
from .reference import Reference


def callable_spans(pos_df: pd.DataFrame, ref: Reference) -> pd.DataFrame:
    """Per-chromosome first/last SNP, SNP count, callable span, and full length."""
    g = pos_df.groupby("chr")["pos0"]
    out = pd.DataFrame({
        "first_snp": g.min(), "last_snp": g.max(), "n_snps": g.size(),
    }).reset_index()
    out["span_bp"] = (out["last_snp"] - out["first_snp"]).clip(lower=0)
    out["full_bp"] = out["chr"].map(ref.core_chrom_lengths_bp).fillna(0).astype(int)
    return out.sort_values("chr")


def per_pair_fraction(blocks_path, sep, bp_per_cm, callable_cm,
                      min_block_snp=IBD_MIN_BLOCK_SNP,
                      min_block_kb=IBD_MIN_BLOCK_KB) -> pd.DataFrame:
    """Per-pair total/max IBD and callable-denominator f. Every pair emitted."""
    df = read_blocks(blocks_path, sep=sep)          # blocks come back half-open
    s1, s2 = df["sample1"].astype(str), df["sample2"].astype(str)
    df["pair"] = np.where(s1 < s2, s1 + "__" + s2, s2 + "__" + s1)

    all_pairs = df[["pair"]].drop_duplicates()

    ibd = df[df["different"] == 0] if "different" in df.columns else df
    # every compared pair is already in `all_pairs`, so filtering here removes spurious
    # short segments from the numerator without losing a pair
    ibd = filter_ibd_blocks(ibd, min_snp=min_block_snp, min_kb=min_block_kb).copy()
    ibd["seg_bp"] = (ibd["end"].astype(int) - ibd["start"].astype(int)).clip(lower=0)
    agg = (ibd.groupby("pair")["seg_bp"]
              .agg(total_ibd_bp="sum", max_ibd_bp="max")
              .reset_index())

    out = (all_pairs.merge(agg, on="pair", how="left")
                    .fillna({"total_ibd_bp": 0, "max_ibd_bp": 0}))
    out["total_ibd_cm"] = out["total_ibd_bp"] / bp_per_cm
    out["max_ibd_cm"] = out["max_ibd_bp"] / bp_per_cm
    out["f"] = out["total_ibd_cm"] / callable_cm
    f_clip = out["f"].clip(lower=1e-6, upper=1.0)
    out["gen_to_mrca_approx"] = np.where(
        out["total_ibd_bp"] > 0, np.log2(1.0 / f_clip).clip(lower=0.0), np.nan)
    return out.sort_values("f", ascending=False).reset_index(drop=True)


def snp_density(pos_df, ref: Reference, bp_per_cm, callable_cm, full_cm, min_snp):
    n_snps = len(pos_df)
    win_bp = bp_per_cm  # 1 cM window
    recs = []
    for chrom, grp in pos_df.groupby("chr"):
        full_bp = ref.core_chrom_lengths_bp.get(chrom)
        if full_bp is None:
            continue
        n_win = int(np.ceil(full_bp / win_bp))
        win_idx = (grp["pos0"].values // win_bp).astype(int)
        counts = np.bincount(win_idx, minlength=n_win)[:n_win]
        for w, c in enumerate(counts):
            recs.append({"chr": chrom, "window_cm": w, "n_snps": int(c)})
    win_df = pd.DataFrame(recs)
    counts = win_df["n_snps"].values
    summary = {
        "n_snps": n_snps,
        "callable_cm": round(callable_cm, 1),
        "full_genome_cm": round(full_cm, 1),
        "snp_per_cm_callable": round(n_snps / callable_cm, 3),
        "snp_per_cm_full_genome": round(n_snps / full_cm, 3),
        "n_windows_1cm": len(win_df),
        "median_snp_per_window": float(np.median(counts)),
        "mean_snp_per_window": round(float(np.mean(counts)), 3),
        "frac_windows_below_min_snp": round(float(np.mean(counts < min_snp)), 4),
        "min_snp_floor": min_snp,
    }
    return summary, win_df
