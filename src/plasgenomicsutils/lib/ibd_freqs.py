"""Global + per-group alternate-allele frequencies from a BCF/VCF.

Alt AF = alt-allele-count / non-missing-allele-count, per ``chr:pos``. The global
table and every group table are computed in a single pass over the file,
accumulating per-group counts as it goes. Group table is group-major with SNP
order following record order.

Genotypes are read as whole per-record numpy arrays (cyvcf2), so allele counting is
a vectorized reduction over all samples rather than a per-sample Python loop.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_allele_freqs(
    bcf_path: str,
    sample_to_group: dict[str, str] | None = None,
    zero_based: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Single pass over ``bcf_path``.

    Parameters
    ----------
    sample_to_group:
        Mapping of sample name -> group for samples present in the metadata.
        Samples absent from the mapping contribute to the global AF only (they
        are excluded from every group). ``None`` computes global AF only.
    zero_based:
        If True, emit ``chr:pos-1`` snp_ids (0-based) to match matrix labels
        from :mod:`plasgenomicsutils.lib.ibd_matrix`.

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
        pos = v.POS - 1 if zero_based else v.POS
        snp_id = f"{v.CHROM}:{pos}"

        # (n_samples, ploidy+1) int; last column is phase, missing allele = -1
        alleles = v.genotype.array()[:, :-1]
        called = alleles >= 0
        is_alt = alleles > 0
        g_an = int(called.sum())
        g_ac = int(is_alt.sum())

        global_rows.append({"snp_id": snp_id, "af": g_ac / g_an if g_an else float("nan")})
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
