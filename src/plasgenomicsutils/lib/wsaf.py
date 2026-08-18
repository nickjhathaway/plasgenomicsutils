"""Within-sample allele frequencies: which polyclonal samples have one dominant clone?

Fws says how clonal a sample is. It does not say whether a sample that fails the gate can
still be *used*, and that is a separate question with a concrete answer: an infection whose
dominant clone carries most of the parasitaemia can be re-genotyped to that clone and treated
as monoclonal, while two strains of comparable size cannot -- a within-sample minor-allele
filter has to stay below 0.5 to keep the dominant call at all, so forcing one there would
delete a real strain and leave a chimera.

So the question this answers is the operational one: **for a sample failing the Fws gate, is
there a dominant clone, and at what fraction?** A dominant clone at 70% means the minor
alleles sit at 30% or below, which is exactly ``filter_ad_regenotype --min-freq 0.30``. The
two numbers are the same choice, so the threshold a caller picks is a statement about
parasite composition rather than a tuning knob.

The evidence is in the read counts. At a heterozygous site the fraction of reads carrying
each allele follows the strain proportions, so for ``K`` strains the het sites fall into bands
at the fractions those proportions imply (Zhu et al. 2019; the basis of DEploid). Filtering at
``f`` zeroes every minor allele below ``f``, so the sites left heterozygous afterwards are
exactly those whose minor fraction is at or above ``f``. Their **rate over the covered sites**
is therefore a direct prediction of what the filter will leave behind, and
``min_freq_needed`` inverts it: the smallest threshold that gets that residue below
``max_residual_het``.

A rate over covered sites is the right denominator, and a share of het sites is not: a sample
with 250 het sites of which 13% are near 0.5 has 33 such sites in a 20,000-site panel, which
is a few repetitive or mismapped loci rather than a second genome. Two genomes in one host
differ across the genome.

Two per-site definitions are reported, because they answer different questions:

``minor_frac``
    the *per-sample* minor fraction: the reads carrying the **second** most-supported allele,
    over the site's depth. In [0, 0.5], anchored to whatever dominates that sample -- which
    is what "is there a dominant clone" is asking, and needs no population. Taking the
    second allele rather than everything-but-the-first keeps the range at 0.5 where a site
    carries three alleles, and makes it the largest single competitor's share.
``wsmaf``
    the within-sample frequency of the **population-level** minor allele, in [0, 1] --
    the quantity in the COI literature (Paschalidis et al. 2023), comparable across
    studies. Needs the population frequencies, so it is computed in the same pass.

Neither is a COI estimate. ``coiaf`` and ``THE REAL McCOIL`` estimate the number of strains;
this predicts whether re-genotyping to the dominant clone will work, and stops there.

References
----------
Zhu, S. J. et al. (2019) The origins and relatedness structure of mixed infections vary
with local prevalence of *P. falciparum* malaria. *eLife* 8, e40845.

Paschalidis, A. et al. (2023) coiaf: directly estimating complexity of infection with
allele frequencies. *PLOS Computational Biology* 19(6), e1010247.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: A sample counts as carrying a minor allele at a site once its fraction clears this.
#: Below it the signal is indistinguishable from sequencing error at ordinary depths.
WSAF_MIN_MINOR = 0.02

#: ...and once the minor allele has at least this many reads, so a single stray read at low
#: depth is not a heterozygous site. Matches filter_ad_regenotype's --min-reads default, so
#: "present" means the same thing in both.
WSAF_MIN_MINOR_READS = 2

#: Share of the parasitaemia the dominant clone must hold for the sample to be worth
#: re-genotyping to it. The complement is the filter threshold: 0.70 here means
#: ``filter_ad_regenotype --min-freq 0.30``. This is a statement about what composition you
#: are willing to call one clone, so it belongs to the caller rather than to this code.
WSAF_MIN_DOMINANT = 0.70

#: Fraction of covered sites allowed to remain heterozygous after that filter before the
#: sample counts as a real mixture rather than a dominant clone with companions.
WSAF_MAX_RESIDUAL_HET = 0.02

#: Covered sites a sample needs before it is judged at all: below this, these rates are
#: counts small enough to turn on a few bad loci.
WSAF_MIN_SITES = 1000

#: Resolution of the per-sample minor-fraction histogram every summary is derived from.
#: Fixed rather than storing each site, so memory is O(samples) and not O(het calls).
WSAF_BIN = 0.01

_EDGES = np.round(np.arange(0.0, 0.5 + WSAF_BIN / 2, WSAF_BIN), 10)
_CENTRES = (_EDGES[:-1] + _EDGES[1:]) / 2
_N_BINS = len(_EDGES) - 1


def _bin_index(frac: float) -> int:
    """Index of the first histogram bin at or above ``frac``."""
    return int(min(_N_BINS, max(0, int(np.ceil(round(frac / WSAF_BIN, 6))))))


def _quantile_of(counts: np.ndarray, q: float) -> float:
    """Quantile of the binned minor fractions (bin centres, so accurate to half a bin)."""
    total = counts.sum()
    if not total:
        return float("nan")
    return float(_CENTRES[int(np.searchsorted(np.cumsum(counts), q * total, side="left"))])


def _mode_of(counts: np.ndarray) -> float:
    """Centre of the fullest bin -- a robust mode for a banded distribution."""
    if not counts.sum():
        return float("nan")
    return float(_CENTRES[int(np.argmax(counts))])


def _count_bands(counts: np.ndarray, min_share: float = 0.10) -> int:
    """How many separated peaks the minor fractions form.

    A band is a run of bins each holding at least ``min_share`` of the fullest bin; runs have
    to be separated by at least one empty-ish bin to count twice. More bands means more
    strain-proportion combinations, i.e. more strains -- though this counts bands, not
    strains, and stops short of inverting one into the other.
    """
    if not counts.sum():
        return 0
    # pairs of bins, so a band is not split in two by one thin bin at this resolution
    coarse = counts[: (_N_BINS // 2) * 2].reshape(-1, 2).sum(axis=1)
    tall = coarse >= max(1.0, coarse.max() * min_share)
    return int(np.diff(np.concatenate(([0], tall.view(np.int8)))).clip(min=0).sum())


def rate_at_or_above(counts: np.ndarray, n_sites: int, frac: float) -> float:
    """Fraction of covered sites whose minor allele reaches ``frac``.

    Equivalently, the share of the genome that stays heterozygous after
    ``filter_ad_regenotype --min-freq frac``, since that filter zeroes everything below it.
    """
    if not n_sites:
        return float("nan")
    return float(counts[_bin_index(frac):].sum() / n_sites)


def min_freq_needed(counts: np.ndarray, n_sites: int,
                    max_residual_het: float = WSAF_MAX_RESIDUAL_HET,
                    min_minor: float = WSAF_MIN_MINOR) -> float:
    """Smallest ``--min-freq`` that leaves under ``max_residual_het`` of sites heterozygous.

    The residue falls as the threshold rises, so this is the first threshold that clears it.
    ``NaN`` when no threshold below 0.5 does -- the sample carries a co-dominant strain, and
    a filter cannot separate the two without deleting one.

    Returns ``min_minor`` when the sample is already under the bar with no filtering at all,
    which is the monoclonal case.
    """
    if not n_sites:
        return float("nan")
    lo = _bin_index(min_minor)
    suffix = np.cumsum(counts[::-1])[::-1] / n_sites
    for k in range(lo, _N_BINS):          # bin k starts at k * WSAF_BIN, the last is 0.49
        if suffix[k] < max_residual_het:
            return float(round(k * WSAF_BIN, 10))
    return float("nan")


def classify_profile(freq_needed: float, n_sites: int,
                     min_sites: int = WSAF_MIN_SITES,
                     min_dominant: float = WSAF_MIN_DOMINANT,
                     min_minor: float = WSAF_MIN_MINOR) -> str:
    """Name a sample by what it would take to reduce it to one clone.

    ``"monoclonal"``
        already under the residual bar with no filtering. Not the same as no heterozygous
        calls: a deeply sequenced clonal sample accumulates a scattering of them from error
        alone, which is why a rate is asked for rather than a count.
    ``"dominant_clone"``
        a filter at or below ``1 - min_dominant`` reduces it to its dominant clone, so it can
        be treated as monoclonal at the composition the caller asked for.
    ``"mixed"``
        it would take a harder filter than that, or none below 0.5 works at all. Either way
        the second strain is too large to remove without deleting it.
    ``"undetermined"``
        too few covered sites to judge. Usually thin coverage.
    """
    if n_sites < min_sites:
        return "undetermined"
    if not np.isfinite(freq_needed):
        return "mixed"
    if freq_needed <= min_minor:
        return "monoclonal"
    if freq_needed <= 1.0 - min_dominant + 1e-9:
        return "dominant_clone"
    return "mixed"


def wsaf_profile(bcf_path: str, *, min_depth: int = 10,
                 min_minor: float = WSAF_MIN_MINOR,
                 min_minor_reads: int = WSAF_MIN_MINOR_READS,
                 min_dominant: float = WSAF_MIN_DOMINANT,
                 max_residual_het: float = WSAF_MAX_RESIDUAL_HET,
                 min_sites: int = WSAF_MIN_SITES,
                 sites_out: str | None = None) -> pd.DataFrame:
    """One row per sample saying whether it has a dominant clone, and at what fraction.

    Single pass over the file. Each sample's minor fractions are accumulated into a fixed
    histogram rather than stored, so memory does not grow with the number of heterozygous
    calls; ``sites_out`` streams to disk as the pass proceeds.

    Parameters
    ----------
    min_depth:
        Per-sample read depth a site needs before its fractions are used. Below this the
        fraction is quantised too coarsely to place a band (at 10x the possible values are
        multiples of 0.1).
    min_minor, min_minor_reads:
        What counts as heterozygous: the minor allele needs both this fraction and this many
        reads. The read floor matters at low depth, where one stray read is 10% of a 10x site.
    min_dominant:
        Share of the parasitaemia the dominant clone must hold. ``1 - min_dominant`` is the
        ``filter_ad_regenotype --min-freq`` this implies.
    max_residual_het:
        Fraction of covered sites allowed to stay heterozygous after that filter.
    min_sites:
        Covered sites needed before a sample is judged rather than left ``"undetermined"``.
    sites_out:
        Optional path for the per-(sample, site) fractions behind the summary, for plotting
        the distributions rather than reading their summaries.

    Returns
    -------
    A tibble-like frame, one row per sample:

    ``n_sites``, ``n_het``, ``het_rate``
        sites passing ``min_depth``, of which how many were heterozygous.
    ``min_freq_needed``
        the smallest ``filter_ad_regenotype --min-freq`` that reduces this sample to one
        clone, or ``NaN`` if none below 0.5 does. **The operational column**: it is the
        argument to pass, per sample.
    ``dominant_frac``
        ``1 - min_freq_needed``, i.e. the share of the parasitaemia the dominant clone holds
        once the rest is removed. Compare it against what you are willing to call one clone.
    ``residual_het_rate``
        fraction of covered sites that would stay heterozygous at the implied threshold, so
        what the filter leaves behind at the composition asked for.
    ``minor_mode``
        the fullest band of the minor fraction, which estimates the minor strain's
        proportion, since that proportion is what sets the band. Where a sample's het calls
        are mostly low-level noise there is no band and the mode lands on the noise floor
        instead, so read it with ``n_bands`` rather than alone.
    ``minor_median``, ``minor_q95``
        quantiles of the same. Descriptive: a quantile is a share of het sites, so on a
        sample with few of them a handful of bad loci move it, which is why the class is
        decided on rates instead.
    ``n_bands``
        separated peaks in the minor fractions; more implies more strains.
    ``wsmaf_mean``
        mean within-sample frequency of the population-level minor allele at het sites,
        the quantity used by COI methods.
    ``class``
        ``monoclonal`` / ``dominant_clone`` / ``mixed`` / ``undetermined``.
    """
    import contextlib

    from cyvcf2 import VCF

    from ..utils.small_utils import Utils

    vcf = VCF(bcf_path)
    samples = list(vcf.samples)
    n = len(samples)

    n_sites = np.zeros(n, dtype=np.int64)
    hist = np.zeros((n, _N_BINS), dtype=np.int64)
    wsmaf_sum = np.zeros(n)
    n_het = np.zeros(n, dtype=np.int64)

    with contextlib.ExitStack() as stack:
        sites_fh = None
        if sites_out:
            sites_fh = stack.enter_context(Utils.smart_open_write(sites_out))
            sites_fh.write("sample\tsnp_id\tminor_frac\talt_frac\tplaf\twsmaf\n")
        stack.callback(vcf.close)
        for v in vcf:
            ad = v.format("AD")
            if ad is None or ad.shape[1] < 2:
                continue
            ad = np.where(ad < 0, 0, np.asarray(ad, dtype=float))
            depth = ad.sum(axis=1)
            ok = depth >= min_depth
            if not ok.any():
                continue
            n_sites += ok

            # the two best-supported alleles per sample; the runner-up's share is the minor
            # fraction, so three alleles at a third each read as 0.33 rather than 0.67
            part = np.partition(ad, -2, axis=1) if ad.shape[1] > 2 else ad
            second = np.sort(part[:, -2:], axis=1)[:, 0]
            mf = np.zeros(n)
            mf[ok] = second[ok] / depth[ok]
            # the population alt frequency this record shows, over the samples with depth --
            # the anchor for wsmaf, computed from the same reads rather than an INFO field
            af_alt = np.zeros(n)
            af_alt[ok] = ad[ok, 1:].sum(axis=1) / depth[ok]
            plaf = float(af_alt[ok].mean())

            het = ok & (mf >= min_minor) & (second >= min_minor_reads)
            if not het.any():
                continue
            idx = np.flatnonzero(het)
            # WSMAF: the population-level minor allele is the alternate where PLAF <= 0.5,
            # otherwise the reference, so the within-sample frequency flips with it
            wsmaf = af_alt[idx] if plaf <= 0.5 else 1.0 - af_alt[idx]
            bins = np.clip((mf[idx] / WSAF_BIN).astype(np.int64), 0, _N_BINS - 1)
            np.add.at(hist, (idx, bins), 1)
            wsmaf_sum[idx] += wsmaf
            n_het[idx] += 1

            if sites_fh is not None:
                snp = f"{v.CHROM}:{v.POS - 1}"
                for j, i in enumerate(idx):
                    sites_fh.write(f"{samples[i]}\t{snp}\t{mf[i]:.6g}\t{af_alt[i]:.6g}\t"
                                   f"{plaf:.6g}\t{wsmaf[j]:.6g}\n")
    rows = []
    for i, s in enumerate(samples):
        counts, cov = hist[i], int(n_sites[i])
        need = min_freq_needed(counts, cov, max_residual_het=max_residual_het,
                              min_minor=min_minor)
        rows.append({
            "sample": s,
            "n_sites": cov,
            "n_het": int(n_het[i]),
            "het_rate": n_het[i] / cov if cov else float("nan"),
            "min_freq_needed": need,
            "dominant_frac": 1.0 - need if np.isfinite(need) else float("nan"),
            "residual_het_rate": rate_at_or_above(counts, cov, 1.0 - min_dominant),
            "minor_mode": _mode_of(counts),
            "minor_median": _quantile_of(counts, 0.5),
            "minor_q95": _quantile_of(counts, 0.95),
            "n_bands": _count_bands(counts),
            "wsmaf_mean": wsmaf_sum[i] / n_het[i] if n_het[i] else float("nan"),
            "class": classify_profile(need, cov, min_sites=min_sites,
                                      min_dominant=min_dominant, min_minor=min_minor),
        })
    out = pd.DataFrame(rows)
    # worst first: the mixtures, then by how hard a filter each one needs
    return (out.sort_values(["residual_het_rate", "min_freq_needed"], ascending=[False, False])
               .reset_index(drop=True))
