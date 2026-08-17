"""Global + per-group alternate-allele frequencies from a BCF/VCF.

Alt AF = alt-allele-count / non-missing-allele-count, per ``chr:pos0``. The global
table and every group table are computed in a single pass over the file,
accumulating per-group counts as it goes. Group table is group-major with SNP
order following record order.

Three readings of "how common is this allele" are reported side by side, because in a
polyclonal infection they genuinely differ:

* ``af``  -- allele count over called alleles. A genotype is a hard call, so a sample
  carrying an allele at 5% within-host counts the same as one carrying it at 95%.
* ``af_weighted`` -- the mean of each sample's *within-sample* frequency (from
  ``FORMAT/AD``), over samples rather than over alleles. Infections share haplotypes, so
  averaging the within-host frequencies tracks the population frequency more closely than
  counting hard calls does.
* ``prevalence`` -- the fraction of samples whose *genotype* carries the allele, which is
  the number usually quoted for a resistance marker.
* ``prevalence_ad`` -- the fraction whose *reads* support it, at ``filter_ad_regenotype``'s
  thresholds. It catches a minor clone the caller left out of the genotype, so it reads
  higher than ``prevalence`` wherever the cohort is polyclonal.

Genotypes are read as whole per-record numpy arrays (cyvcf2), so allele counting is
a vectorized reduction over all samples rather than a per-sample Python loop.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .intervals import snp_label



def _within_sample_fractions(v):
    """((fractions, depths), missing) -- AD / sum(AD) per sample, and the raw AD counts.

    The denominator is that sample's own read depth, so each sample contributes a
    frequency between 0 and 1 whatever its coverage -- a 500x sample does not outvote a
    50x one, which is the point of averaging over samples rather than over reads.
    """
    try:
        ad = v.format("AD")
    except KeyError:
        return None, True
    if ad is None:
        return None, True
    ad = np.asarray(ad, dtype=float)
    if ad.ndim != 2 or ad.shape[1] != len(v.ALT) + 1:
        return None, True
    ad[ad < 0] = 0.0                      # cyvcf2 writes missing as a large negative
    tot = ad.sum(axis=1)
    out = np.full(ad.shape, np.nan)
    ok = tot > 0
    out[ok] = ad[ok] / tot[ok, None]
    return (out, ad), not ok.any()


def _freq_row(hit, called, sample_called, frac, k, mask, ad_min_reads=2,
              ad_min_freq=0.01):
    """One frequency row for the samples in `mask`.

    `hit` marks the allele slots that are the allele of interest; `k` is its index in the
    record's alleles (None when every ALT is being taken together).

    Two prevalences, because "carries this allele" has two readings. `prevalence` counts
    samples whose *genotype* has it -- what the caller committed to. `prevalence_ad` counts
    samples whose *reads* support it at `ad_min_reads` and `ad_min_freq`, which catches a
    minor clone the caller did not call. The thresholds are `filter_ad_regenotype`'s, so
    "present" means the same thing here as it does there.
    """
    an = int(called[mask].sum())
    ac = int(hit[mask].sum())
    af = ac / an if an else float("nan")

    n_samples = int(sample_called[mask].sum())
    n_alt = int((hit[mask] & called[mask]).any(axis=1).sum())

    af_w = float("nan")
    n_ad = 0
    n_alt_ad = 0
    prev_ad = float("nan")
    if frac is not None:
        f, depth = frac[0][mask], frac[1][mask]
        # every ALT together is one minus the reference fraction, which keeps a
        # multiallelic site consistent with its per-ALT rows summing to the same total
        col = 1.0 - f[:, 0] if k is None else f[:, k]
        reads = depth[:, 1:].sum(axis=1) if k is None else depth[:, k]
        usable = ~np.isnan(col)
        n_ad = int(usable.sum())
        if n_ad:
            af_w = float(col[usable].sum() / n_ad)
            present = usable & (reads >= ad_min_reads) & (col >= ad_min_freq)
            n_alt_ad = int(present.sum())
            prev_ad = n_alt_ad / n_ad

    return {
        "af": af, "maf": min(af, 1 - af), "ac": ac, "an": an,
        "af_weighted": af_w, "n_samples_ad": n_ad,
        "prevalence": n_alt / n_samples if n_samples else float("nan"),
        "n_samples_alt": n_alt, "n_samples": n_samples,
        "prevalence_ad": prev_ad, "n_samples_alt_ad": n_alt_ad,
    }


def compute_allele_freqs(
    bcf_path: str,
    sample_to_group: dict[str, str] | None = None,
    with_pos_vcf: bool = False,
    per_alt: bool = False,
    weighted: bool = True,
    ad_min_reads: int = 2,
    ad_min_freq: float = 0.01,
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
    per_alt:
        One row per (SNP, ALT allele) instead of one per SNP, with ``alt`` and
        ``alt_index`` naming which allele the row is about. Off by default, since a
        table with one row per SNP is what the selection statistic joins on. At a
        biallelic site the two forms agree; at a multiallelic one the default collapses
        every ALT together and this separates them.
    weighted:
        Compute ``af_weighted`` and ``prevalence_ad`` from ``FORMAT/AD`` (default). Costs a
        second array read per record; turn it off for speed, or where AD is absent -- in
        which case both are NaN and a note says so once.
    ad_min_reads, ad_min_freq:
        What counts as the allele being present in a sample for ``prevalence_ad``: at least
        this many reads *and* this within-sample fraction (defaults 2 and 0.01, matching
        ``filter_ad_regenotype``, so "present" means the same thing in both).

    ``snp_id`` is always the canonical 0-based ``chr:pos0`` label
    (:func:`~plasgenomicsutils.lib.intervals.snp_label`), matching the IBD matrix
    columns; the record's own ``ID`` field is ignored.

    Returns
    -------
    (global_df, group_df). Columns: ``snp_id`` (and ``group`` on the group table),
    ``af``, ``maf``, ``ac``, ``an``, ``af_weighted``, ``n_samples_ad``, ``prevalence``,
    ``n_samples_alt``, ``n_samples``, ``prevalence_ad``, ``n_samples_alt_ad``; plus
    ``alt``/``alt_index`` when ``per_alt``.
    ``ac``/``an`` are allele counts, ``n_samples_*`` are sample counts -- at ploidy 2
    they differ by a factor of two, which is why both are reported.
    group_df is empty if no mapping is given.
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

    warned_no_ad = False

    for v in vcf:
        pos0 = v.POS - 1                        # VCF is 1-based; everything inward is not
        snp_id = snp_label(v.CHROM, pos0)

        # (n_samples, ploidy+1) int; last column is phase, missing allele = -1
        alleles = v.genotype.array()[:, :-1]
        called = alleles >= 0
        sample_called = called.any(axis=1)

        # within-sample allele fractions from AD, aligned to the record's alleles
        frac = None
        if weighted:
            frac, missing_ad = _within_sample_fractions(v)
            if missing_ad and not warned_no_ad:
                print(f"  note: no usable FORMAT/AD at {snp_id} (and possibly elsewhere); "
                      f"af_weighted is NaN there")
                warned_no_ad = True

        # which allele each row is about: every ALT together, or one row per ALT
        if per_alt:
            targets = [(k, v.ALT[k - 1]) for k in range(1, len(v.ALT) + 1)]
        else:
            targets = [(None, None)]

        for k, alt in targets:
            hit = alleles == k if k is not None else alleles > 0
            extra = {} if k is None else {"alt": alt, "alt_index": k}
            grow = _freq_row(hit, called, sample_called, frac, k,
                             np.ones(len(samples), bool), ad_min_reads, ad_min_freq)
            grow = {"snp_id": snp_id, **extra, **grow}
            if with_pos_vcf:
                grow["pos_vcf"] = v.POS
            global_rows.append(grow)
            for r in groups:
                m = group_masks[r]
                row = _freq_row(hit, called, sample_called, frac, k, m,
                                ad_min_reads, ad_min_freq)
                group_rows[r].append({"group": r, "snp_id": snp_id, **extra, **row})

    vcf.close()

    global_df = pd.DataFrame(global_rows)
    if groups:
        group_df = pd.concat([pd.DataFrame(group_rows[r]) for r in groups], ignore_index=True)
    else:
        cols = ["group", "snp_id"] + (["alt", "alt_index"] if per_alt else []) + [
            "af", "maf", "ac", "an", "af_weighted", "n_samples_ad",
            "prevalence", "n_samples_alt", "n_samples",
            "prevalence_ad", "n_samples_alt_ad"]
        group_df = pd.DataFrame(columns=cols)
    return global_df, group_df
