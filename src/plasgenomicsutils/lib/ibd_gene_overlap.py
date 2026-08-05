"""Per-gene IBD block overlap between sample groups.

For each gene interval and each pair of groups, the fraction of sample pairs that share
an IBD block **overlapping the gene** -- not merely pairs that are IBD at a genotyped SNP
inside it. A pair counts when any of its IBD segments overlaps the (optionally padded)
gene interval, so a long segment spanning a gene is captured even with no panel SNP inside.

The denominator is the total number of pairs compared in that group-pair, so pairs that
are never IBD still count against the fraction: the analyzed-sample set is taken from
**every** row of the blocks file (hmmibd-rs emits a segment for every compared pair), then
only ``different == 0`` (IBD) segments drive the numerator.
"""

from __future__ import annotations

from itertools import combinations_with_replacement

import numpy as np
import pandas as pd

from .reference import normalise_chr


def _analyzed_group_counts(blocks_df: pd.DataFrame, sample_to_group: dict) -> dict:
    """Count analyzed samples per group (samples appearing anywhere in the blocks)."""
    analyzed = set(blocks_df["sample1"]).union(blocks_df["sample2"])
    counts: dict = {}
    for s in analyzed:
        g = sample_to_group.get(s)
        if g is not None and not (isinstance(g, float) and np.isnan(g)):
            counts[g] = counts.get(g, 0) + 1
    return counts


def gene_block_overlap(
    blocks_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    sample_to_group: dict,
    within: int = 0,
    only_ibd: bool = True,
) -> pd.DataFrame:
    """Per-gene, per-group-pair IBD-block overlap fraction.

    blocks_df : hmmibd-rs blocks (``sample1, sample2, chr, start, end[, different]``).
    genes_df  : gene intervals (``name, chr, start, end``).
    sample_to_group : sample -> group label.
    within    : pad each gene interval by this many bp on both sides (default 0).

    Coordinates are 0-based half-open ``[start, end)`` (BED). ``hmmibd-rs`` reports both
    block endpoints as the 0-based position of the segment's first and last SNP, so pass
    blocks through :func:`blocks_to_half_open` first to make ``end`` exclusive.
    Returns a data frame:
    ``gene, chr, start, end, group_a, group_b, n_pairs_ibd, n_pairs_total, frac_pairs_ibd``
    (one row per gene x group-pair; every group-pair is emitted, even with 0 sharing).
    """
    counts = _analyzed_group_counts(blocks_df, sample_to_group)
    groups = sorted(counts)
    if len(groups) < 2:
        raise SystemExit("need at least two groups among the analyzed samples")

    def total_pairs(a: str, b: str) -> int:
        return counts[a] * (counts[a] - 1) // 2 if a == b else counts[a] * counts[b]

    ibd = blocks_df
    if only_ibd and "different" in blocks_df.columns:
        ibd = blocks_df[blocks_df["different"] == 0]
    ibd = ibd.loc[:, ["sample1", "sample2", "chr", "start", "end"]].copy()
    ibd["chr"] = ibd["chr"].map(normalise_chr)
    by_chr = {c: d.reset_index(drop=True) for c, d in ibd.groupby("chr", sort=False)}

    has_gid = "gene_id" in genes_df.columns          # GFF Names repeat; gene_id is unique
    gp_index = list(combinations_with_replacement(groups, 2))   # (a, b) with a <= b
    rows = []
    for _, gene in genes_df.iterrows():
        chrom = normalise_chr(gene["chr"])
        lo = int(gene["start"]) - within
        hi = int(gene["end"]) + within
        gp_counts = {gp: 0 for gp in gp_index}
        d = by_chr.get(chrom)
        if d is not None and len(d):
            m = (d["start"].to_numpy() < hi) & (d["end"].to_numpy() > lo)
            if m.any():
                s1 = d["sample1"].to_numpy()[m]
                s2 = d["sample2"].to_numpy()[m]
                seen: set = set()
                for a_s, b_s in zip(s1, s2):
                    pair = (a_s, b_s) if a_s < b_s else (b_s, a_s)
                    if pair in seen:
                        continue
                    seen.add(pair)
                    ga = sample_to_group.get(a_s)
                    gb = sample_to_group.get(b_s)
                    if ga is None or gb is None or ga not in counts or gb not in counts:
                        continue
                    key = (ga, gb) if ga <= gb else (gb, ga)
                    gp_counts[key] += 1
        for a, b in gp_index:
            tp = total_pairs(a, b)
            num = gp_counts[(a, b)]
            row = {"gene": gene["name"]}
            if has_gid:
                row["gene_id"] = gene["gene_id"]
            row.update({
                "chr": chrom, "start": int(gene["start"]), "end": int(gene["end"]),
                "group_a": a, "group_b": b,
                "n_pairs_ibd": num, "n_pairs_total": tp,
                "frac_pairs_ibd": (num / tp) if tp else float("nan"),
            })
            rows.append(row)
    return pd.DataFrame(rows)
