"""Build / load the binary (pairs x SNPs) IBD matrix from hmmibd-rs blocks.

The matrix is assembled from COO triplets and a single ``.tocsr()``. Overlapping
blocks for a pair collapse back to a binary 1 (``.data[:] = 1``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, load_npz, save_npz

from .vcf_io import SnpPanel


BLOCKS_DTYPE = {
    "sample1": str, "sample2": str, "chr": str, "start": int, "end": int,
}


def read_blocks(path: str, sep: str = "\t") -> pd.DataFrame:
    """Read an hmmibd-rs blocks TSV (sample1, sample2, chr, start, end, different, Nsnp)."""
    return pd.read_csv(path, sep=sep, comment="#", dtype=BLOCKS_DTYPE)


def build_pair_index(blocks_df: pd.DataFrame) -> tuple[dict, list]:
    """Enumerate order-normalised unique pairs -> (pair_to_row, pair_labels)."""
    pairs = set()
    for s1, s2 in zip(blocks_df["sample1"], blocks_df["sample2"]):
        pairs.add((s1, s2) if s1 < s2 else (s2, s1))
    pair_list = sorted(pairs)
    pair_to_row = {p: i for i, p in enumerate(pair_list)}
    pair_labels = [f"{p[0]}__{p[1]}" for p in pair_list]
    return pair_to_row, pair_labels


def build_matrix(
    blocks_df: pd.DataFrame,
    panel: SnpPanel,
    pair_to_row: dict,
    only_different_zero: bool = True,
) -> csr_matrix:
    """Fill a binary (n_pairs x n_snps) CSR matrix from IBD blocks.

    only_different_zero: keep only ``different == 0`` blocks (true IBD);
    pass ``False`` for ``--all-blocks``.
    """
    n_pairs = len(pair_to_row)
    n_snps = len(panel)

    if only_different_zero and "different" in blocks_df.columns:
        work_df = blocks_df[blocks_df["different"] == 0]
    else:
        work_df = blocks_df

    rows_acc: list[np.ndarray] = []
    cols_acc: list[np.ndarray] = []
    for row in work_df.itertuples(index=False):
        s1, s2 = row.sample1, row.sample2
        key = (s1, s2) if s1 < s2 else (s2, s1)
        pair_row = pair_to_row.get(key)
        if pair_row is None:
            continue
        col_indices = panel.snps_in_block(row.chr, int(row.start), int(row.end))
        if col_indices.size:
            cols_acc.append(col_indices)
            rows_acc.append(np.full(col_indices.size, pair_row, dtype=np.int64))

    if rows_acc:
        rows = np.concatenate(rows_acc)
        cols = np.concatenate(cols_acc)
        data = np.ones(rows.size, dtype=np.uint8)
        csr = coo_matrix((data, (rows, cols)), shape=(n_pairs, n_snps), dtype=np.uint8).tocsr()
        csr.data[:] = 1  # collapse summed duplicates from overlapping blocks -> binary
    else:
        csr = csr_matrix((n_pairs, n_snps), dtype=np.uint8)
    return csr


def save_matrix(mat: csr_matrix, pair_labels: list, snp_labels: list, out_prefix: str) -> None:
    save_npz(f"{out_prefix}.npz", mat)
    Path(f"{out_prefix}.pair_labels.txt").write_text("\n".join(pair_labels) + "\n")
    Path(f"{out_prefix}.snp_labels.txt").write_text("\n".join(snp_labels) + "\n")


def load_matrix(out_prefix: str) -> tuple[csr_matrix, list, list]:
    # accept the prefix (as build_ibd_matrix writes) or the full "<prefix>.npz" path
    if out_prefix.endswith(".npz"):
        out_prefix = out_prefix[: -len(".npz")]
    mat = load_npz(f"{out_prefix}.npz")
    pairs = Path(f"{out_prefix}.pair_labels.txt").read_text().splitlines()
    snps = Path(f"{out_prefix}.snp_labels.txt").read_text().splitlines()
    return mat, pairs, snps


def pair_summary(mat: csr_matrix, pair_labels: list) -> pd.DataFrame:
    """Per-pair total IBD SNPs and fraction of the SNP panel in IBD."""
    n_snps = mat.shape[1]
    counts = np.asarray(mat.sum(axis=1)).ravel()
    return pd.DataFrame({
        "pair": pair_labels,
        "n_ibd_snps": counts,
        "frac_ibd": counts / n_snps,
    })
