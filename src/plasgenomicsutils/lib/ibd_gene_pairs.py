"""Sample pairs sharing IBD over each gene, with how much of the gene they cover.

The pair-level companion to :mod:`.ibd_gene_overlap`: instead of a fraction per
group-pair, one row per sample pair x IBD block x gene. Pairs with no IBD over a gene are
simply absent -- use the overlap table when the denominator (all pairs compared) matters.

A pair appears more than once for a gene only if it has several separate IBD segments
spanning it. Intervals are 0-based half-open (see :mod:`.intervals`).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .ibd_matrix import IBD_MIN_BLOCK_KB, IBD_MIN_BLOCK_SNP, filter_ibd_blocks
from .intervals import overlaps
from .reference import normalise_chr


COLUMNS = [
    "sample1", "sample2", "chr", "block_start", "block_end",
    "gene", "gene_id", "gene_start", "gene_end",
    "coverage", "covered_start", "covered_end", "covered_bp", "percent_covered",
    "gene_cluster_id", "gene_cluster_size",
]


def single_linkage(a, b):
    """Connected components over an edge list, as ``(id, size)`` maps keyed by sample.

    Single-linkage clustering: a sample joins a cluster if it shares IBD with **any**
    member, so a chain of pairs is one cluster even where its ends never share directly.
    Ids run largest cluster first, so ``1`` is the biggest group at that gene and the
    numbering means the same thing from one gene to the next. Ties break on the first
    member's name, which keeps the ids stable between runs.
    """
    nodes = sorted(set(a) | set(b))
    index = {n: i for i, n in enumerate(nodes)}
    parent = list(range(len(nodes)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]      # path compression
            i = parent[i]
        return i

    for x, y in zip(a, b):
        rx, ry = find(index[x]), find(index[y])
        if rx != ry:
            parent[rx] = ry

    members: dict[int, list[str]] = {}
    for n in nodes:
        members.setdefault(find(index[n]), []).append(n)
    order = sorted(members, key=lambda r: (-len(members[r]), members[r][0]))
    cid = {r: i + 1 for i, r in enumerate(order)}
    ids = {n: cid[find(index[n])] for n in nodes}
    sizes = {n: len(members[find(index[n])]) for n in nodes}
    return ids, sizes


def gene_ibd_pairs(
    blocks_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    *,
    within: int = 0,
    only_ibd: bool = True,
    min_block_snp: int = IBD_MIN_BLOCK_SNP,
    min_block_kb: float = IBD_MIN_BLOCK_KB,
) -> pd.DataFrame:
    """One row per sample pair x IBD block x gene.

    blocks_df : hmmibd-rs blocks (``sample1, sample2, chr, start, end[, different]``),
        already half-open (see :func:`.intervals.blocks_to_half_open`, applied by
        :func:`.ibd_matrix.read_blocks`).
    genes_df  : gene intervals (``name, chr, start, end``; ``gene_id`` optional).
    within    : pad each gene when deciding overlap. Coverage is always measured against
        the gene's own span, so a block reaching only into the padding covers 0.
    """
    ibd = blocks_df
    if only_ibd and "different" in blocks_df.columns:
        ibd = blocks_df[blocks_df["different"] == 0]
    ibd = filter_ibd_blocks(ibd, min_snp=min_block_snp, min_kb=min_block_kb)
    ibd = ibd.loc[:, ["sample1", "sample2", "chr", "start", "end"]].copy()
    ibd["chr"] = ibd["chr"].map(normalise_chr)
    by_chr = {c: d.reset_index(drop=True) for c, d in ibd.groupby("chr", sort=False)}

    has_gid = "gene_id" in genes_df.columns
    frames = []
    for _, gene in genes_df.iterrows():
        chrom = normalise_chr(gene["chr"])
        gs, ge = int(gene["start"]), int(gene["end"])
        d = by_chr.get(chrom)
        if d is None or not len(d):
            continue
        m = overlaps(d["start"].to_numpy(), d["end"].to_numpy(), gs - within, ge + within)
        if not m.any():
            continue
        d = d[m]

        bs = d["start"].to_numpy(); be = d["end"].to_numpy()
        cs = np.maximum(bs, gs)                      # intersection with the gene itself
        ce = np.minimum(be, ge)
        covered = np.maximum(0, ce - cs)
        complete = covered >= (ge - gs)
        s1 = np.minimum(d["sample1"].to_numpy(), d["sample2"].to_numpy())
        s2 = np.maximum(d["sample1"].to_numpy(), d["sample2"].to_numpy())

        frames.append(pd.DataFrame({
            "sample1": s1, "sample2": s2, "chr": chrom,
            "block_start": bs, "block_end": be,
            "gene": str(gene["name"]),
            "gene_id": str(gene["gene_id"]) if has_gid else pd.NA,
            "gene_start": gs, "gene_end": ge,
            "coverage": np.where(complete, "complete", "partial"),
            "covered_start": np.where(covered > 0, cs, np.nan),
            "covered_end": np.where(covered > 0, ce, np.nan),
            "covered_bp": covered,
            "percent_covered": 100.0 * covered / (ge - gs) if ge > gs else np.nan,
        }))
        ids, sizes = single_linkage(s1, s2)
        frames[-1]["gene_cluster_id"] = [ids[v] for v in s1]
        frames[-1]["gene_cluster_size"] = [sizes[v] for v in s1]

    if not frames:
        return pd.DataFrame(columns=COLUMNS)
    out = pd.concat(frames, ignore_index=True)
    out = out.sort_values(["gene", "sample1", "sample2"]).reset_index(drop=True)
    return out.loc[:, COLUMNS]
