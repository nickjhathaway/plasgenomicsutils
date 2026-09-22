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
    counters = {"n_variants": 0, "n_low_depth": 0, "n_star_alleles": 0, "n_star_only": 0}

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
            # Per ALT allele, not per record. At a triallelic site where one sample carries
            # ALT1 and another ALT2 there are two singletons; counting carriers of *any*
            # alternate saw two non-reference samples and booked a doubleton shared between
            # them -- which then fed the "near-identical to X" flag with evidence of
            # near-identity between samples carrying different alleles.
            alts = list(v.ALT)
            n_alt = len(alts)
            if not n_alt:
                continue
            # `*` is a spanning deletion, not an alternate base, and a singleton count is
            # about alleles. Counting it does two kinds of damage: a sample that is the only
            # one with a deletion is booked a private SNP it does not have, and -- the
            # expensive one -- two samples that share a deletion are booked a DOUBLETON,
            # which feeds `shared[(i, j)]` and therefore the "near-identical to X" flag.
            # Deletions are haplotype markers, so sharing one is common and the fabricated
            # evidence accumulates between samples that are not near-identical at all.
            real = [i for i, a in enumerate(alts, start=1) if a != "*"]
            if not real:
                counters["n_star_only"] += 1
                continue
            counters["n_star_alleles"] += n_alt - len(real)
            if n_alt == 1:
                # the fast path, and bit-for-bit what this always did on a biallelic record
                per_allele = [((gt == 1) | (gt == 3)) & ~miss]
            else:
                alleles = v.genotype.array()[:, :-1]      # (samples, ploidy) allele indices
                per_allele = [np.isin(alleles, a).any(axis=1) & ~miss for a in real]
            for carriers in per_allele:
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
    df.attrs["n_star_alleles"] = counters["n_star_alleles"]
    df.attrs["n_star_only_records"] = counters["n_star_only"]
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


#: Counts returned by :func:`singleton_to_missing`.
SINGLETON_RECODE_COUNTS = (
    "records", "records_with_singleton_alt", "records_all_alts_singleton",
    "alts_recoded", "records_with_singleton_star", "stars_recoded",
    "calls_recoded", "calls_fully_missing",
)

#: What :func:`singleton_to_missing` is allowed to recode.
SINGLETON_RECODE_MODES = ("real", "star", "all")


def singleton_allele_indices(alts, allele_array, *, min_samples: int = 1,
                             real: bool = True, star: bool = True) -> list[int]:
    """Which ALT indices (1-based) are singletons that this mode may blank.

    The one place the rule lives, so a caller re-genotyping in memory
    (:func:`~plasgenomicsutils.lib.regenotype.filter_ad_regenotype`) and one rewriting a file
    (:func:`singleton_to_missing`) cannot drift apart.

    ``alts`` is the ALT list and ``allele_array`` the ``(samples, ploidy)`` array of allele
    indices, missing slots negative. A sample carries allele *k* if its genotype names *k*
    anywhere, so a ``1/2`` het is one carrier of each.

    ``star`` is honoured only where a real alternate survives to be saved by the ``*``
    going: one nobody is blanking and somebody carries. Blanking a ``*`` earns its keep by
    rescuing a SNP from ``--snps-only``; with no such SNP there is nothing to rescue, the
    deletion is simply what the record holds, and removing it would empty the ALT column for
    nothing. See :data:`~plasgenomicsutils.lib.allele_counts.ALT_SAMPLE_MAX_TAG`, which makes
    the same carve-out for the carrier count.
    """
    real_idx = [i + 1 for i, a in enumerate(alts) if a != "*"]
    star_idx = [i + 1 for i, a in enumerate(alts) if a == "*"]
    if not (real_idx if real else []) and not (star_idx if star else []):
        return []
    carriers = {k: int(np.isin(allele_array, k).any(axis=1).sum())
                for k in real_idx + star_idx}
    out = [k for k in real_idx if carriers[k] <= min_samples] if real else []
    if star and any(k not in out and carriers[k] > 0 for k in real_idx):
        out += [k for k in star_idx if carriers[k] <= min_samples]
    return out


def singleton_to_missing(inp: str, out: str, *, min_samples: int = 1,
                         alleles: str = "real") -> dict[str, int]:
    """Null every genotype slot naming an ALT carried by ``<= min_samples`` samples.

    The per-allele counterpart of the record-level singleton drop, and the singleton
    analogue of :func:`~plasgenomicsutils.lib.spanning_del.spanning_del_to_missing`. A
    record-level filter keeps the whole record when *some* alternate has enough carriers, so
    a singleton alternate sitting beside a well-supported one rides through: the record is
    kept and the singleton allele stays on it. This blanks the calls that name such an
    alternate instead, so a following ``--trim-alt-alleles`` removes the allele and the
    record keeps only its supported ones.

    The ALT column is left exactly as it was, so this is separable from the trim and can be
    checked on its own -- same contract as ``spanning_del_to_missing``.

    ``alleles`` says which alternates are eligible: ``"real"`` (the default) the bases and
    indels only, ``"star"`` the spanning deletion ``*`` only, ``"all"`` both.

    ``*`` is kept apart because it is not a variant allele: it says a deletion called
    elsewhere covers the position (which is why ``spanning_del_filter`` exists and is
    default-off), so it is never counted as *support* for a record. But a ``*`` carried by one
    sample is still a private observation, and left on a record it takes the site out of any
    SNP panel built with ``--snps-only`` -- a common SNP lost for one sample's deletion. The
    ``"star"`` mode is the singleton-only counterpart of ``spanning_del_to_missing``: it
    blanks the ``*`` calls of that one sample and nothing else, so the deletion is never
    treated as support and a *common* deletion is left exactly as it was, for
    ``spanning_del_filter`` to judge.

    A ``*`` is recoded **only where a real alternate survives to be saved** -- one that is
    not itself being recoded and that some sample carries. Blanking a ``*`` is worth doing
    because it rescues a SNP from ``--snps-only``; where there is no such SNP there is nothing
    to rescue, and the deletion is simply what the record holds. That covers the record whose
    sole ALT is ``*`` (the same carve-out
    :data:`~plasgenomicsutils.lib.allele_counts.ALT_SAMPLE_MAX_TAG` makes for the carrier
    count: there the deletion is the variant and its own count decides it) and the record
    whose real alternates are left uncarried. Either way the record is written through with
    its ``*`` intact, so the recode can never be the thing that empties an ALT column --
    which matters where something puts a dropped record back, since a whitelist rescue would
    otherwise restore an allele-less record rather than the deletion it was asked to keep.

    A carrier count is per sample: a sample carries allele
    *k* if its genotype names *k* anywhere, so a ``1/2`` het is one carrier of each -- the
    same rule :func:`count_singletons` and
    :func:`~plasgenomicsutils.lib.allele_counts.add_alt_sample_counts` use.

    Counts returned: ``records`` seen; ``records_with_singleton_alt`` (at least one real ALT
    recoded); ``records_all_alts_singleton`` (every real ALT was a singleton, so the record
    becomes ref-only once trimmed -- judged over the real alternates whatever ``alleles``
    says, since a surviving ``*`` is not a variant); ``alts_recoded`` (real ALT alleles nulled
    across all records); ``records_with_singleton_star`` and ``stars_recoded`` (the same for
    ``*``); ``calls_recoded`` (sample-calls that lost at least one slot); and
    ``calls_fully_missing`` (calls left with nothing).
    """
    if alleles not in SINGLETON_RECODE_MODES:
        raise ValueError(f"alleles must be one of {SINGLETON_RECODE_MODES}, got {alleles!r}")
    do_real = alleles in ("real", "all")
    do_star = alleles in ("star", "all")
    from cyvcf2 import VCF, Writer

    vcf = VCF(inp)
    writer = Writer(out, vcf)
    st = dict.fromkeys(SINGLETON_RECODE_COUNTS, 0)
    try:
        for v in vcf:
            st["records"] += 1
            alts = list(v.ALT)
            real = [i + 1 for i, a in enumerate(alts) if a != "*"]
            arr = v.genotype.array()[:, :-1]     # (samples, ploidy) allele indices
            singleton = singleton_allele_indices(alts, arr, min_samples=min_samples,
                                                 real=do_real, star=do_star)
            if not singleton:
                writer.write_record(v)
                continue

            real_hit = [k for k in singleton if k in real]
            star_hit = [k for k in singleton if k not in real]
            if real_hit:
                st["records_with_singleton_alt"] += 1
                st["alts_recoded"] += len(real_hit)
                if len(real_hit) == len(real):
                    st["records_all_alts_singleton"] += 1
            if star_hit:
                st["records_with_singleton_star"] += 1
                st["stars_recoded"] += len(star_hit)

            drop = set(singleton)
            gts = v.genotypes                    # [[a1, a2, ..., phased], ...]
            changed = False
            for g in gts:
                hit = [k for k in range(len(g) - 1) if g[k] in drop]
                if not hit:
                    continue
                for k in hit:
                    g[k] = -1
                st["calls_recoded"] += 1
                if all(a < 0 for a in g[:-1]):
                    st["calls_fully_missing"] += 1
                changed = True
            if changed:
                v.genotypes = gts
            writer.write_record(v)
    finally:
        writer.close()
        vcf.close()
    return st
