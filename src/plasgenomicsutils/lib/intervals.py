"""Genomic coordinate conventions and helpers.

**Everything in this package is 0-based.** Intervals are half-open ``[start, end)`` -- the
BED convention -- so ``end - start`` is the width in bp and intervals that merely touch
(``end1 == start2``) do not overlap. Variant positions are 0-based too, and a variant at
``pos0`` lies inside ``[start, end)`` when ``start <= pos0 < end``. There is one rule; no
part of the package is 1-based.

Formats that use something else are converted once, at the boundary, and their numbering
never propagates inward:

* **VCF/BCF** ``POS`` is 1-based, as the format defines: it becomes ``pos0 = POS - 1`` the
  moment it is read.
* **hmmibd-rs** blocks report ``start`` and ``end`` as the 0-based positions of the first
  and last SNP of a segment, both inclusive. :func:`blocks_to_half_open` shifts ``end`` by
  one so block intervals compose with everything else.

SNP identity follows from that: :func:`snp_label` is the *only* way a ``snp_id`` is built,
always from ``(chrom, pos0)``. Ids found in an input file -- a BED name column, a VCF
``ID`` field set by ``bcftools annotate --set-id`` -- are never adopted as keys, because
whether such an id used ``%POS`` or ``%POS0`` is unknowable from the file. They are carried
alongside as ``source_id`` for traceability only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# Stamped into outputs that carry SNP labels so a consumer can verify what it was given
# rather than infer it from the numbers.
SNP_COORD_SYSTEM = "0-based"


def check_snp_coord_system(found, source: str) -> None:
    """Fail loudly when a file's SNP labels are not in the current coordinate system.

    Labels alone cannot be trusted to reveal a mismatch: shift a dense panel by one and
    most labels still land on a real neighbouring SNP, so the join succeeds and silently
    pairs each SNP with the wrong allele frequency. The stamp is what makes it detectable.
    """
    if found == SNP_COORD_SYSTEM:
        return
    what = "carries no coordinate-system stamp" if found is None else f"is '{found}'"
    raise SystemExit(
        f"ERROR: {source} {what}, but this version writes and expects "
        f"'{SNP_COORD_SYSTEM}' SNP ids (chr:pos0).\n"
        "  Files written before the switch to a single 0-based convention used 1-based "
        "ids. Regenerate them -- do not mix, because an off-by-one panel still joins "
        "against neighbouring SNPs and would corrupt the result without erroring."
    )


def snp_label(chrom, pos0):
    """The canonical ``chr:pos0`` SNP id. The only place a ``snp_id`` is constructed."""
    if isinstance(chrom, (str, bytes)) or np.isscalar(chrom):
        return f"{chrom}:{int(pos0)}"
    return [f"{c}:{int(p)}" for c, p in zip(chrom, pos0)]


def vcf_pos(pos0):
    """The 1-based VCF position for a 0-based coordinate (for human cross-reference)."""
    return np.asarray(pos0) + 1


def blocks_to_half_open(blocks_df: pd.DataFrame, *, end_col: str = "end") -> pd.DataFrame:
    """Return a copy of hmmibd-rs blocks with ``end`` made exclusive.

    hmmibd-rs emits the 0-based position of the segment's last SNP as ``end``; adding one
    turns the segment into a half-open ``[start, end)`` interval.
    """
    out = blocks_df.copy()
    out[end_col] = out[end_col].astype("int64") + 1
    return out


def overlaps(
    starts: np.ndarray, ends: np.ndarray, lo: int, hi: int
) -> np.ndarray:
    """Boolean mask of half-open intervals ``[starts, ends)`` overlapping ``[lo, hi)``."""
    return (np.asarray(starts) < hi) & (np.asarray(ends) > lo)


def overlap_span(start: int, end: int, lo: int, hi: int) -> tuple[int, int]:
    """Intersection of two half-open intervals; empty when ``start >= end``."""
    return max(int(start), int(lo)), min(int(end), int(hi))


def variant_in_interval(pos0: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """Mask of 0-based variant positions falling inside the half-open ``[lo, hi)``."""
    p = np.asarray(pos0)
    return (p >= lo) & (p < hi)
