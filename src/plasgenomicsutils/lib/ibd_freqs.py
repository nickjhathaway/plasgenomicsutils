"""Global + per-group alternate-allele frequencies from a BCF/VCF.

Alt AF = alt-allele-count / non-missing-allele-count, per ``chr:pos0``. The global
table and every group table are computed in a single pass over the file,
accumulating per-group counts as it goes. Group table is group-major with SNP
order following record order.

Genotypes are read as whole per-record numpy arrays (cyvcf2), so allele counting is
a vectorized reduction over all samples rather than a per-sample Python loop.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .intervals import snp_label


def compute_allele_freqs(
    bcf_path: str,
    sample_to_group: dict[str, str] | None = None,
    with_pos_vcf: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Single pass over ``bcf_path``.

    Parameters
    ----------
    sample_to_group:
        Mapping of sample name -> group for samples present in the metadata.
        Samples absent from the mapping contribute to the global AF only (they
        are excluded from every group). ``None`` computes global AF only.
    with_pos_vcf:
        Add a ``pos_vcf`` column holding the 1-based VCF position, for looking a
        variant up by eye. Off by default -- it is derivable from ``snp_id`` and
        only inflates the file.

    ``snp_id`` is always the canonical 0-based ``chr:pos0`` label
    (:func:`~plasgenomicsutils.lib.intervals.snp_label`), matching the IBD matrix
    columns; the record's own ``ID`` field is ignored.

    Returns
    -------
    (global_df, group_df) where global_df has columns [snp_id, af] and group_df
    has columns [group, snp_id, af]. group_df is empty if no mapping is given.
    """
    from cyvcf2 import VCF

    vcf = VCF(bcf_path)
    samples = list(vcf.samples)

    groups: list[str] = []
    group_of_sample: dict[str, str] = {}
    if sample_to_group:
        group_of_sample = {s: sample_to_group[s] for s in samples if s in sample_to_group}
        groups = sorted(set(group_of_sample.values()))
    # boolean sample masks (aligned to the file's sample order), one per group
    group_masks = {
        r: np.fromiter((group_of_sample.get(s) == r for s in samples), dtype=bool, count=len(samples))
        for r in groups
    }

    global_rows: list[dict] = []
    # group -> list of {group, snp_id, af} rows, kept in record order
    group_rows: dict[str, list[dict]] = {r: [] for r in groups}

    for v in vcf:
        pos0 = v.POS - 1                        # VCF is 1-based; everything inward is not
        snp_id = snp_label(v.CHROM, pos0)

        # (n_samples, ploidy+1) int; last column is phase, missing allele = -1
        alleles = v.genotype.array()[:, :-1]
        called = alleles >= 0
        is_alt = alleles > 0
        g_an = int(called.sum())
        g_ac = int(is_alt.sum())

        grow = {"snp_id": snp_id, "af": g_ac / g_an if g_an else float("nan")}
        if with_pos_vcf:
            grow["pos_vcf"] = v.POS
        global_rows.append(grow)
        for r in groups:
            m = group_masks[r]
            an = int(called[m].sum())
            ac = int(is_alt[m].sum())
            group_rows[r].append({
                "group": r, "snp_id": snp_id,
                "af": ac / an if an else float("nan"),
            })

    vcf.close()

    global_df = pd.DataFrame(global_rows)
    if groups:
        group_df = pd.concat([pd.DataFrame(group_rows[r]) for r in groups], ignore_index=True)
    else:
        group_df = pd.DataFrame(columns=["group", "snp_id", "af"])
    return global_df, group_df
