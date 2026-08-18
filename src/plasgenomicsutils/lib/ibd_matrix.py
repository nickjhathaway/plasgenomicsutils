"""Build / load the binary (pairs x SNPs) IBD matrix from hmmibd-rs blocks.

The matrix is assembled from COO triplets and a single ``.tocsr()``. Overlapping
blocks for a pair collapse back to a binary 1 (``.data[:] = 1``).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, load_npz, save_npz

from .intervals import SNP_COORD_SYSTEM, blocks_to_half_open, check_snp_coord_system
from .vcf_io import SnpPanel


BLOCKS_DTYPE = {
    "sample1": str, "sample2": str, "chr": str, "start": int, "end": int,
}


def read_blocks(path: str, sep: str = "\t") -> pd.DataFrame:
    """Read an hmmibd-rs blocks TSV (sample1, sample2, chr, start, end, different, Nsnp).

    hmmibd-rs reports ``start`` and ``end`` as the 0-based positions of the segment's first
    and last SNP, both inclusive; ``end`` is shifted by one so the returned intervals are
    half-open ``[start, end)`` like every other interval here (see :mod:`.intervals`).
    """
    df = pd.read_csv(path, sep=sep, comment="#", dtype=BLOCKS_DTYPE)
    return blocks_to_half_open(df)


# Short, SNP-poor IBD segments are commonly spurious, so they are dropped by default
# rather than left for every caller to pre-filter. These are the conventional thresholds.
IBD_MIN_BLOCK_SNP = 15
IBD_MIN_BLOCK_KB = 15.0


def filter_ibd_blocks(
    blocks_df: pd.DataFrame,
    min_snp: int = IBD_MIN_BLOCK_SNP,
    min_kb: float = IBD_MIN_BLOCK_KB,
) -> pd.DataFrame:
    """Drop IBD segments shorter than ``min_kb`` or carrying fewer than ``min_snp`` SNPs.

    Apply this to the **IBD evidence only** (``different == 0`` rows). The set of pairs that
    were compared -- the denominator in every fraction -- must still come from the
    unfiltered frame, or dropping a pair's only short segment would also drop the pair.

    ``0`` (or ``None``) disables either criterion. Blocks are half-open, so a segment's
    length is ``end - start``.
    """
    if not len(blocks_df):
        return blocks_df
    keep = pd.Series(True, index=blocks_df.index)
    if min_kb:
        span = blocks_df["end"].astype("int64") - blocks_df["start"].astype("int64")
        keep &= span >= float(min_kb) * 1000.0
    if min_snp:
        if "Nsnp" not in blocks_df.columns:
            warnings.warn(
                f"--min-block-snp {min_snp} was requested but the blocks have no 'Nsnp' "
                "column; only the length filter was applied",
                stacklevel=2,
            )
        else:
            keep &= blocks_df["Nsnp"].astype("int64") >= int(min_snp)
    return blocks_df[keep]


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
    min_block_snp: int = IBD_MIN_BLOCK_SNP,
    min_block_kb: float = IBD_MIN_BLOCK_KB,
) -> csr_matrix:
    """Fill a binary (n_pairs x n_snps) CSR matrix from IBD blocks.

    only_different_zero: keep only ``different == 0`` blocks (true IBD);
    pass ``False`` for ``--all-blocks``.

    Short / SNP-poor segments are dropped first (see :func:`filter_ibd_blocks`); the pair
    index, built by the caller from the unfiltered frame, still carries every compared
    pair, so the denominator is unaffected.
    """
    n_pairs = len(pair_to_row)
    n_snps = len(panel)

    if only_different_zero and "different" in blocks_df.columns:
        work_df = blocks_df[blocks_df["different"] == 0]
    else:
        work_df = blocks_df
    work_df = filter_ibd_blocks(work_df, min_snp=min_block_snp, min_kb=min_block_kb)

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
    """Write the sparse pairs x SNPs matrix plus its row and column labels to ``path``."""
    save_npz(f"{out_prefix}.npz", mat)
    Path(f"{out_prefix}.pair_labels.txt").write_text("\n".join(pair_labels) + "\n")
    # The coordinate system is recorded rather than left to be inferred from the numbers:
    # a panel whose labels are off by one still *matches* neighbouring SNPs in a dense
    # panel, so a mismatch cannot be detected reliably by comparing label sets.
    Path(f"{out_prefix}.snp_labels.txt").write_text(
        f"#snp_coord_system={SNP_COORD_SYSTEM}\n" + "\n".join(snp_labels) + "\n")


def load_matrix(out_prefix: str) -> tuple[csr_matrix, list, list]:
    """Read back what :func:`save_matrix` wrote: ``(matrix, pair_labels, snp_labels)``."""
    # accept the prefix (as build_ibd_matrix writes) or the full "<prefix>.npz" path
    if out_prefix.endswith(".npz"):
        out_prefix = out_prefix[: -len(".npz")]
    mat = load_npz(f"{out_prefix}.npz")
    pairs = Path(f"{out_prefix}.pair_labels.txt").read_text().splitlines()
    lines = Path(f"{out_prefix}.snp_labels.txt").read_text().splitlines()
    check_snp_coord_system(_read_coord_header(lines), f"{out_prefix}.snp_labels.txt")
    snps = [ln for ln in lines if not ln.startswith("#")]
    return mat, pairs, snps


def _read_coord_header(lines: list) -> str | None:
    for ln in lines[:1]:
        if ln.startswith("#snp_coord_system="):
            return ln.split("=", 1)[1].strip()
    return None


def pair_summary(mat: csr_matrix, pair_labels: list) -> pd.DataFrame:
    """Per-pair total IBD SNPs and fraction of the SNP panel in IBD."""
    n_snps = mat.shape[1]
    counts = np.asarray(mat.sum(axis=1)).ravel()
    return pd.DataFrame({
        "pair": pair_labels,
        "n_ibd_snps": counts,
        "frac_ibd": counts / n_snps,
    })
