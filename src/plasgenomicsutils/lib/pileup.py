"""Shared single-site pileup primitives for read-level diagnostics.

Both strand-bias tools ask the same low-level question — for one position, which
reads cover it, on which strand, with which base at what sequencing cycle — so
that iteration and the true-cycle calculation live here once.

The **true sequencing cycle** matters: reverse reads are stored reverse-complemented
relative to the reference, so a base at query index ``q`` in a length-``L`` reverse
read was actually sequenced at cycle ``L-1-q``. Miscalls clustered at one cycle
implicate a flow-cell/phasing defect; spread across the read body implicates
strand-specific chemistry.
"""

from __future__ import annotations

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Reverse-complement a nucleotide string (N preserved)."""
    return seq.translate(_COMP)[::-1]


def parse_pos(spec: str):
    """Split a ``CHROM:POS`` (1-based) string into ``(chrom, pos1)``."""
    chrom, pos = spec.rsplit(":", 1)
    return chrom, int(pos)


def lookup_ref_base(ref_fasta: str, chrom: str, pos0: int) -> str:
    """Fetch the single uppercase reference base at 0-based ``pos0``."""
    import pysam

    fa = pysam.FastaFile(ref_fasta)
    try:
        return fa.fetch(chrom, pos0, pos0 + 1).upper()
    finally:
        fa.close()


def iter_site_reads(bam, chrom: str, pos0: int, *, min_mapq: int = 0):
    """Yield ``(pileup_read, alignment, base, seq_cycle)`` for reads over ``pos0``.

    ``pos0`` is 0-based. Secondary/supplementary/unmapped/duplicate/qcfail reads,
    reads failing ``min_mapq``, and reads with no aligned base at the site
    (deletion/ref-skip) are skipped — the same primary-read tally both tools use.
    ``base`` is uppercase; ``seq_cycle`` is the true sequencing cycle.
    """
    for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, stepper="nofilter",
                          min_base_quality=0, max_depth=10_000_000,
                          ignore_overlaps=False):
        if col.reference_pos != pos0:
            continue
        for pr in col.pileups:
            read = pr.alignment
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate or read.is_qcfail:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            seq = read.query_sequence
            if seq is None:
                continue
            q = pr.query_position
            base = seq[q].upper()
            length = len(seq)
            seq_cycle = (length - 1 - q) if read.is_reverse else q
            yield pr, read, base, seq_cycle
