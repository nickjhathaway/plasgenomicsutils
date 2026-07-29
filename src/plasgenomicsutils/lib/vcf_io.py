"""Shared VCF/BCF + BED readers used across the IBD and filtering layers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .reference import normalise_chr
from ..utils.small_utils import Utils


class SnpPanel:
    """A SNP panel loaded from a VCF or BED, plus a fast per-chromosome index.

    Column convention:
      * ``snp_id``  — ``"chr:pos"`` label (VCF uses 1-based POS; BED uses start)
      * ``chr``     — chromosome string as written in the source file
      * ``pos0``    — 0-based coordinate (VCF POS-1; BED start)

    The 0-based ``pos0`` axis is what hmmibd-rs blocks (0-based inclusive) index
    into, so the IBD matrix columns and these labels line up.
    """

    def __init__(self, df: pd.DataFrame):
        self.df = df.reset_index(drop=True)
        self.labels: list[str] = self.df["snp_id"].tolist()
        self._index = self._build_index(self.df)

    # -- loaders ------------------------------------------------------------
    @classmethod
    def from_vcf(cls, path: str) -> "SnpPanel":
        """Load SNPs from a VCF/BCF text stream (1-based POS -> 0-based pos0)."""
        rows = []
        with Utils.smart_open_read(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t", 5)
                chrom, pos1 = parts[0], int(parts[1])
                rows.append({"snp_id": f"{chrom}:{pos1}", "chr": chrom, "pos0": pos1 - 1})
        return cls(pd.DataFrame(rows))

    @classmethod
    def from_bed(cls, path: str) -> "SnpPanel":
        """Load SNPs from a BED (0-based start used as the SNP coordinate)."""
        rows = []
        with Utils.smart_open_read(path) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                parts = line.rstrip("\n").split("\t")
                chrom, start = parts[0], int(parts[1])
                snp_id = parts[3] if len(parts) >= 4 and parts[3] else f"{chrom}:{start}"
                rows.append({"snp_id": snp_id, "chr": chrom, "pos0": start})
        return cls(pd.DataFrame(rows))

    @classmethod
    def load(cls, path: str, fmt: str) -> "SnpPanel":
        if fmt == "vcf":
            return cls.from_vcf(path)
        if fmt == "bed":
            return cls.from_bed(path)
        raise SystemExit(f"ERROR: unknown --snp-format '{fmt}' (expected vcf|bed)")

    # -- index / lookup -----------------------------------------------------
    @staticmethod
    def _build_index(snp_df: pd.DataFrame) -> dict:
        """chr -> {"pos0": sorted np.array, "global_idx": column-position array}."""
        index = {}
        for chrom, grp in snp_df.groupby("chr"):
            grp_sorted = grp.sort_values("pos0")
            index[chrom] = {
                "pos0": grp_sorted["pos0"].values,
                "global_idx": grp_sorted.index.values,
            }
        return index

    def snps_in_block(self, chrom: str, start0: int, end0_inclusive: int) -> np.ndarray:
        """Global SNP column indices with 0-based pos in ``[start0, end0]`` on ``chrom``.

        Binary search -> O(log n + k). ``chrom`` must match the source spelling.
        """
        if chrom not in self._index:
            return np.array([], dtype=np.int64)
        pos = self._index[chrom]["pos0"]
        gidx = self._index[chrom]["global_idx"]
        lo = np.searchsorted(pos, start0, side="left")
        hi = np.searchsorted(pos, end0_inclusive, side="right")
        return gidx[lo:hi]

    def __len__(self) -> int:
        return len(self.df)


def positions_frame(path: str, fmt: str) -> pd.DataFrame:
    """Load just (normalised-chr, pos) rows for callable-span / density work.

    The chromosome is normalised (``Pf3D7_07_v3`` -> ``7``) so it keys against
    the reference registry's chromosome lengths.
    """
    rows = []
    with Utils.smart_open_read(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            rows.append((normalise_chr(p[0]), int(p[1])))  # vcf POS or bed start
    return pd.DataFrame(rows, columns=["chr", "pos"])
