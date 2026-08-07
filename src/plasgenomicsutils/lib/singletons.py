"""Per-sample singleton counts -- a genetic-outlier flag for sample QC.

A *singleton* here is a variant where exactly one sample in the cohort carries the
non-reference allele. Every cohort has some, spread thinly across samples, and a sample
that departs from that in **either** direction is worth looking at:

  * **Too many.** Contamination, a mixed or off-species infection, or an alignment
    problem generating private false positives. MalariaGEN drop samples on exactly this
    criterion when assembling the Pf analysis sets.
  * **Too few.** Another copy of the sample is absorbing them. If two samples carry
    near-identical parasites then neither one's private variants are private any more --
    each becomes a *doubleton* shared with the other, and both singleton counts collapse.
    This is easy to miss, because nothing about the sample looks low-quality: its call
    rate is normal, it just has no variation of its own.

So the same scan also records, for every sample, which other sample it shares its
doubletons with, turning "this count looks odd" into a named pair.

**What that pair means is not decided here.** A near-identical pair is either the same
parasite sequenced twice or two hosts infected by the same clone, and nothing in an
allele count can tell those apart -- in a high-transmission setting the second is common
and is a finding rather than an error. The flag says *near-identical*; use IBD, the
collection metadata, or the sequencing records to decide which. The threshold is
deliberately strict (see ``DEFAULT_DUPLICATE_FRAC``): merely *related* samples sit well
below it, and calling them duplicates would throw away real data.

The count on its own is not comparable between samples with different amounts of missing
data, so the rate over called genotypes is reported alongside it and the outlier call is
made on the rate. Outliers are flagged by median absolute deviation rather than mean/SD,
because the outliers being looked for would otherwise inflate the very spread they are
measured against.
"""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

#: Flag a sample above this many MADs from the cohort median singleton rate.
DEFAULT_MAD_CUTOFF = 5.0
#: Share of a sample's doubletons pointing at one partner before calling the two
#: near-identical. Set high on purpose: on a real 249-sample cohort the median share is
#: ~0.04, genuinely near-identical pairs sit above 0.95, and merely closely-related pairs
#: reach ~0.6-0.8 -- which is where a looser threshold starts discarding real samples.
DEFAULT_DUPLICATE_FRAC = 0.9
#: Scale factor making the MAD a consistent estimator of the SD for normal data.
_MAD_TO_SD = 1.4826


def _read_depth(v, n):
    """Per-sample depth for one record: summed AD if present, else DP.

    AD is preferred because a GATK ``DP`` can be large while the depth actually
    supporting the alleles is not, which is the whole reason the AD-based re-genotyping
    step exists. Missing entries come back negative from cyvcf2 and are read as 0.
    """
    ad = v.format("AD")
    if ad is not None:
        return np.clip(np.asarray(ad, dtype=np.int64), 0, None).sum(axis=1)
    dp = v.format("DP")
    if dp is None:
        return np.full(n, -1, dtype=np.int64)      # no depth to judge by
    return np.clip(np.asarray(dp, dtype=np.int64).ravel(), 0, None)


def count_singletons(vcf_path, samples=None, regions=None, max_missing_frac=1.0,
                     min_depth=0):
    """Count, per sample, the variants where it is the only non-reference carrier.

    Args:
        vcf_path: VCF/BCF path (indexed if ``regions`` is used).
        samples: Restrict to these samples. Singleton status is judged **within the
            set analysed**, so subsetting changes what counts as private.
        regions: Iterable of ``chrom`` or ``chrom:start-end`` strings.
        max_missing_frac: Skip variants whose missing-genotype fraction exceeds this;
            a variant called in only a handful of samples makes a private allele
            trivially likely.
        min_depth: Treat a genotype backed by fewer than this many reads as uncalled.
            Callsets built from gVCFs routinely emit ``0/0`` at sites with **no reads
            at all** rather than ``./.``, which would otherwise be counted as a
            confident reference call and inflate every denominator here. Depth is read
            from ``AD`` (summed) when present, else ``DP``.

    Returns:
        ``(DataFrame, n_variants)``. One row per sample: ``sample, n_called,
        n_singleton, singleton_rate`` (singletons per 1000 called genotypes).
    """
    from cyvcf2 import VCF

    vcf = VCF(vcf_path, samples=samples) if samples else VCF(vcf_path)
    names = list(vcf.samples)
    n = len(names)
    if n == 0:
        raise SystemExit(f"{vcf_path}: no samples")
    singles = np.zeros(n, dtype=np.int64)
    doubles = np.zeros(n, dtype=np.int64)
    called = np.zeros(n, dtype=np.int64)
    shared = Counter()                      # (i, j) -> doubletons carried by both
    counters = {"n_variants": 0, "n_low_depth": 0}

    def scan(it):
        for v in it:
            # gt_types: 0 hom-ref, 1 het, 2 unknown, 3 hom-alt (cyvcf2's coding)
            gt = v.gt_types
            miss = gt == 2
            if min_depth > 0:
                depth = _read_depth(v, n)
                shallow = (depth >= 0) & (depth < min_depth) & ~miss
                counters["n_low_depth"] += int(shallow.sum())
                miss = miss | shallow
            if miss.mean() > max_missing_frac:
                continue
            counters["n_variants"] += 1
            np.add(called, ~miss, out=called)
            carriers = ((gt == 1) | (gt == 3)) & ~miss
            k = carriers.sum()
            if k == 1:
                singles[np.argmax(carriers)] += 1
            elif k == 2:
                i, j = np.flatnonzero(carriers)
                doubles[i] += 1
                doubles[j] += 1
                shared[(int(i), int(j))] += 1

    if regions:
        for r in regions:
            scan(vcf(r))
    else:
        scan(vcf)
    vcf.close()

    # the partner each sample shares most of its doubletons with
    best = np.zeros(n, dtype=np.int64)
    best_of = np.full(n, -1, dtype=np.int64)
    for (i, j), c in shared.items():
        if c > best[i]:
            best[i], best_of[i] = c, j
        if c > best[j]:
            best[j], best_of[j] = c, i

    rate = np.divide(1000.0 * singles, called, out=np.zeros(n), where=called > 0)
    frac = np.divide(best.astype(float), doubles, out=np.zeros(n, dtype=float),
                     where=doubles > 0)
    df = pd.DataFrame({
        "sample": names, "n_called": called, "n_singleton": singles,
        "singleton_rate": rate, "n_doubleton": doubles,
        "top_partner": [names[k] if k >= 0 else "" for k in best_of],
        "n_shared_with_partner": best, "frac_doubletons_with_partner": frac})
    df.attrs["n_low_depth"] = counters["n_low_depth"]
    return df, counters["n_variants"]


def flag_outliers(df, mad_cutoff=DEFAULT_MAD_CUTOFF, column="singleton_rate",
                  duplicate_frac=DEFAULT_DUPLICATE_FRAC):
    """Add ``mad_score``, ``outlier`` and ``flag`` columns.

    The flag is two-sided, because the tails mean different things: an excess of private
    variants points at contamination or mis-alignment, a deficit means another sample is
    absorbing them. The near-identical call is made on the doubleton share alone, not on
    the MAD score, so it still works in a cohort with enough duplicated samples to drag
    the median rate down.
    """
    v = df[column].to_numpy(dtype=float)
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med))) * _MAD_TO_SD
    score = np.zeros_like(v) if mad == 0 else (v - med) / mad

    partner = df["top_partner"] if "top_partner" in df else pd.Series([""] * len(df))
    frac = (df["frac_doubletons_with_partner"] if "frac_doubletons_with_partner" in df
            else pd.Series(np.zeros(len(df))))
    flag = []
    for s, p, f in zip(score, partner, frac):
        parts = []
        if s > mad_cutoff:
            parts.append("excess private variants")
        elif s < -mad_cutoff:
            parts.append("few private variants")
        if f >= duplicate_frac and p:
            parts.append(f"near-identical to {p}")
        flag.append("; ".join(parts))
    return df.assign(mad_score=score, outlier=np.abs(score) > mad_cutoff,
                     flag=pd.Series(flag, index=df.index))
