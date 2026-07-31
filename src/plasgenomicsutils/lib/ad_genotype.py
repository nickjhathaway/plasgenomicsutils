"""Shared AD cleaning + re-genotyping for haploid (blood-stage Plasmodium) calls.

Two problems use this core:
  * within-sample cleaning of low-level AD artifacts, then re-genotyping, and
  * harmonizing separately-called cohorts before merging.

Re-genotyping ranks alleles by depth and only calls a heterozygote when the
minor allele's within-sample frequency is >= ``het_min_af``; otherwise the
sample is called homozygous for the major allele. In a haploid organism a real
``0/1`` indicates a mixed (polyclonal) infection, so this threshold suppresses
sequencing-error hets while keeping genuine mixed-infection signal.
"""

from __future__ import annotations

import numpy as np


def regenotype_from_ad(ad_counts, het_min_af: float = 0.2, called_alleles=None):
    """Re-call a genotype from per-allele AD counts.

    Returns a ``(a, b)`` allele-index tuple, or ``None`` when the sample has no
    supporting reads (missing).

      * total 0                    -> None (missing)
      * only one allele supported  -> homozygous for it
      * minor-allele AF < het_min_af -> homozygous for the major allele
      * otherwise                  -> heterozygous for the top two alleles

    ``called_alleles`` restricts the candidates to the allele indices the caller
    already used for this sample. This only narrows an existing call (a hom-ref
    sample can lose support and go missing, but never gain an alternate allele),
    which preserves an upstream caller's genotypes when its likelihood-based
    calls are trusted over raw read counts. Leave it ``None`` to re-derive the
    genotype from AD alone.
    """
    total = sum(ad_counts)
    if total == 0:
        return None

    candidates = list(range(len(ad_counts)))
    if called_alleles is not None:
        allowed = {a for a in called_alleles if a is not None and a >= 0}
        candidates = [i for i in candidates if i in allowed]
        if not candidates:
            return None
        supported = [i for i in candidates if ad_counts[i] > 0]
        if not supported:
            return None
        if len(supported) == 1:
            return (supported[0], supported[0])
        a, b = sorted(supported[:2])
        return (a, b)

    ranked = sorted(((i, ad_counts[i]) for i in candidates), key=lambda x: -x[1])
    top_idx, _ = ranked[0]

    if len(ranked) < 2 or ranked[1][1] == 0:
        return (top_idx, top_idx)

    second_idx, second_count = ranked[1]
    if second_count / total < het_min_af:
        return (top_idx, top_idx)
    a, b = sorted([top_idx, second_idx])
    return (a, b)


def regenotype_matrix(ad, het_min_af: float = 0.2):
    """Vectorized :func:`regenotype_from_ad` over a cleaned (n_samples, n_alleles) AD matrix.

    Returns ``(gt_a, gt_b)`` int arrays of allele indices, with ``-1, -1`` for
    samples with no supporting reads. Reproduces the scalar rule exactly, including
    its tie-breaking (rank by depth; ties by ascending allele index via a stable
    sort). Only the default path — pass ``called_alleles`` sample-by-sample to
    :func:`regenotype_from_ad` for the restrict-to-called-alleles behavior.
    """
    ad = np.asarray(ad, dtype=float)
    m, k = ad.shape
    total = ad.sum(axis=1)
    order = np.argsort(-ad, axis=1, kind="stable")  # ties -> ascending index
    rows = np.arange(m)
    top = order[:, 0]
    if k >= 2:
        second = order[:, 1]
        second_count = ad[rows, second]
    else:
        second = top.copy()
        second_count = np.zeros(m)
    with np.errstate(invalid="ignore", divide="ignore"):
        minor_af = np.where(total > 0, second_count / total, 0.0)
    only_one = (k < 2) | (second_count == 0)
    het = (~only_one) & (minor_af >= het_min_af)
    gt_a = np.where(het, np.minimum(top, second), top)
    gt_b = np.where(het, np.maximum(top, second), top)
    missing = total == 0
    gt_a = np.where(missing, -1, gt_a)
    gt_b = np.where(missing, -1, gt_b)
    return gt_a.astype(int), gt_b.astype(int)


def clean_ad_matrix(ad, depth, min_reads: int, min_freq: float,
                    protect_ref: bool = False):
    """Zero out per-sample allele depths that fall below the thresholds.

    Parameters
    ----------
    ad : (n_samples, n_alleles) array of allele depths.
    depth : (n_samples,) per-sample total depth used as the frequency denominator.
    min_reads : zero an allele whose read count is < this.
    min_freq : zero an allele whose within-sample frequency (count/depth) is < this.
    protect_ref : if True, never zero the reference allele (column 0).

    Returns a new float array with sub-threshold entries set to 0.
    """
    ad = ad.astype(float).copy()
    depth = depth.astype(float)
    start = 1 if protect_ref else 0
    for a in range(start, ad.shape[1]):
        fail_depth = ad[:, a] < min_reads
        with np.errstate(invalid="ignore", divide="ignore"):
            freq = ad[:, a] / depth
        fail_freq = freq < min_freq
        ad[np.where(fail_depth | fail_freq), a] = 0
    return ad
