"""Fws within-host diversity statistic (Manske 2012; the moimix::getFws estimator).

Fws measures how much of a sample's genetic diversity is captured within the host
versus the population — a monoclonal infection has Fws ~ 1, a polyclonal one less.
It needs a cohort (to estimate population allele frequencies) and per-sample allele
depths (AD). With ``k`` alleles at site ``i`` (REF plus every ALT):

  * population allele frequencies:   P_ik = Σ_s ad_sik / Σ_s Σ_k ad_sik
  * population heterozygosity:       Hs_i = 1 - Σ_k P_ik²
  * within-sample allele fractions:  q_sik = ad_sik / Σ_k ad_sik
  * within-sample heterozygosity:    Hw_si = 1 - Σ_k q_sik²
  * sites are binned by the population **minor-allele fraction** — the share of reads
    not on the major allele, ``1 - max_k P_ik`` — into ``n_bins`` equal bins over
    [0, 0.5]; per bin the mean Hw (per sample) and mean Hs are taken.

At a biallelic site these are exactly moimix's quantities (``1 - (p² + q²)`` and
``min(p, 1 - p)``); the multiallelic form is the same heterozygosity written for
``k`` alleles, so a sample carrying two different non-reference alleles is seen as
mixed. Splitting multiallelics first (``bcftools norm -m-``) hides that: each split
record keeps only REF and one ALT, so a 50:50 mix of two ALTs looks homozygous twice.
Readers therefore take unsplit input and collapse multiallelic records by default;
``multiallelic="skip"`` drops them instead (what moimix does on split input).

Two estimators combine the binned means into Fws — see :func:`compute_fws`:

  * ``"regression"`` (default): Fws = 1 - β, where β is the slope of a regression of
    the binned sample-het means on the binned population-het means, forced through
    the origin. This matches ``moimix::getFws``.
  * ``"ratio"``: Fws = 1 - Σ_bins mean(Hw) / Σ_bins mean(Hs) — a simpler ratio of the
    summed binned means.

The two agree in spirit but not to the digit (the regression weights bins by the
squared population het), so pick deliberately and don't mix a threshold tuned on one
with values from the other.

Because heterozygosity is written for ``k`` alleles, a **microhaplotype** locus from an
amplicon panel -- one record whose alleles are the haplotypes and whose depths are the
per-haplotype read counts -- is just another multiallelic site. Such loci often have no
allele above 50%, so ``1 - max(p)`` runs past the biallelic grid's top of 0.5; the grid
extends upward in the same bin width to hold them, so they bin with sites of similar
heterozygosity rather than being clipped in with the 0.45-0.5 biallelic ones. For a panel
that is mostly such loci ``n_bins=0`` regresses per locus instead of per bin (see
:func:`compute_fws`).

Depths are read from a bcftools-query **AD table** (:func:`read_ad_table`), a **VCF/BCF**
(:func:`read_ad_vcf`) or a long-format **allele table** (:func:`read_allele_table`) into an
:class:`AlleleDepths`, which feeds :func:`compute_fws`.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass

import numpy as np

from .reporting import detail, say

#: How a record with more than one ALT is treated by the readers.
MULTIALLELIC_MODES = ("collapse", "skip")


def _check_multiallelic(mode):
    if mode not in MULTIALLELIC_MODES:
        raise ValueError(f"multiallelic must be one of {MULTIALLELIC_MODES}, not {mode!r}")


# --------------------------------------------------------------------------- #
#  Region exclusion (CNV windows depress Fws, so they are dropped)            #
# --------------------------------------------------------------------------- #


def load_exclude_regions(path):
    """Read a ``chrom, call_start, call_end`` TSV into ``{chrom: [(start, end), ...]}``.

    A CNV inside its own call window inflates within-sample heterozygosity and would
    make a monoclonal sample look polyclonal, so those windows are excluded from Fws.
    Returns an empty dict when ``path`` is falsy.
    """
    excl: dict[str, list[tuple[int, int]]] = {}
    if not path:
        return excl
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            excl.setdefault(row["chrom"], []).append(
                (int(row["call_start"]), int(row["call_end"])))
    return excl


def _in_excluded(chrom, pos, excl):
    return any(s <= pos <= e for s, e in excl.get(chrom, ()))


# --------------------------------------------------------------------------- #
#  Allele depths, reduced to what the estimator needs                         #
# --------------------------------------------------------------------------- #


@dataclass
class AlleleDepths:
    """Per-site, per-sample allele depths reduced to what Fws needs.

    Keeping the full ``[n_sites, n_samples, n_alleles]`` array would grow with the most
    multiallelic site in the file; the estimator only ever needs, per site and sample,
    the total depth and the sum of squared allele depths (for ``Hw = 1 - Σ q²``), plus
    the per-site population allele totals (for ``Hs`` and the minor-allele fraction).
    Those are fixed-size whatever ``k`` is.

    ``ref``, ``depth``, ``sumsq`` are float ``[n_sites, n_samples]``; ``pop`` is float
    ``[n_sites, max_alleles]``, zero-padded past a site's own allele count. A biallelic
    site reduces to the classic ref/alt form exactly (``sumsq = ref² + alt²``).
    """

    ref: np.ndarray
    depth: np.ndarray
    sumsq: np.ndarray
    pop: np.ndarray
    n_multiallelic: int = 0
    multiallelic: str = "collapse"
    n_alt_trimmed: int = 0    # records that lost >= 1 ALT for having no reads in the cohort
    n_monomorphic: int = 0    # records dropped because no ALT had any reads in the cohort
    sites: list | None = None      # one id per site ("CHROM:POS", or the locus name)
    alleles: list | None = None    # per site, the allele strings in ``pop`` column order

    @property
    def n_sites(self) -> int:
        return int(self.depth.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.depth.shape[1])

    @property
    def nonref(self) -> np.ndarray:
        """Reads on any non-reference allele: at a biallelic site, the alt depth."""
        return self.depth - self.ref

    @classmethod
    def from_ref_alt(cls, ref, alt) -> "AlleleDepths":
        """Build from classic biallelic ``[n_sites, n_samples]`` ref/alt depth matrices."""
        ref = np.asarray(ref, dtype=float)
        alt = np.asarray(alt, dtype=float)
        if ref.size == 0:
            n = ref.shape[1] if ref.ndim == 2 else 0
            return cls.empty(n)
        pop = np.stack([ref.sum(axis=1), alt.sum(axis=1)], axis=1)
        return cls(ref=ref, depth=ref + alt, sumsq=ref * ref + alt * alt, pop=pop)

    @classmethod
    def from_records(cls, records, n_samples, **counts) -> "AlleleDepths":
        """Build from an iterable of per-site ``[n_samples, n_alleles]`` AD arrays.

        Allele 0 is REF; a site may have any number of alleles. Missing depths must
        already be zero. ``counts`` are the bookkeeping fields (``n_multiallelic``,
        ``multiallelic``, ``n_alt_trimmed``, ``n_monomorphic``) plus, optionally,
        ``sites`` and ``alleles`` -- the ids and allele names that let population
        frequencies from elsewhere be matched up (:func:`compute_fws` ``pop_freqs``) and
        this cohort's be written out (:meth:`population_freqs`).
        """
        ref_rows, depth_rows, sq_rows, pop_rows = [], [], [], []
        for ad in records:
            ad = np.asarray(ad, dtype=float)
            ref_rows.append(ad[:, 0])
            depth_rows.append(ad.sum(axis=1))
            sq_rows.append((ad * ad).sum(axis=1))
            pop_rows.append(ad.sum(axis=0))
        if not ref_rows:
            return cls.empty(n_samples, **counts)
        max_k = max(len(r) for r in pop_rows)
        pop = np.zeros((len(pop_rows), max_k), dtype=float)
        for i, r in enumerate(pop_rows):
            pop[i, :len(r)] = r
        return cls(ref=np.array(ref_rows, dtype=float),
                   depth=np.array(depth_rows, dtype=float),
                   sumsq=np.array(sq_rows, dtype=float), pop=pop, **counts)

    @classmethod
    def empty(cls, n_samples, **counts) -> "AlleleDepths":
        z = np.empty((0, n_samples), dtype=float)
        return cls(ref=z, depth=z.copy(), sumsq=z.copy(), pop=np.empty((0, 2), dtype=float),
                   **counts)

    def population_freqs(self):
        """This cohort's population allele frequencies, as ``{site: {allele: freq}}``.

        Needs ``sites`` and ``alleles``; frequencies are pooled read fractions, the same
        ones :func:`compute_fws` uses when no ``pop_freqs`` are supplied.
        """
        if self.sites is None or self.alleles is None:
            raise ValueError("these depths carry no site ids / allele names")
        out = {}
        for i, (site, names) in enumerate(zip(self.sites, self.alleles)):
            tot = self.pop[i, :len(names)].sum()
            if tot > 0:
                out[site] = {a: float(self.pop[i, k] / tot) for k, a in enumerate(names)}
        return out

    def multiallelic_note(self) -> str:
        """One line on what happened to multiallelic and read-less alleles, or '' if nothing did.

        The multiallelic count is of records with reads on more than one ALT *in this
        cohort*; an ALT that is only there because the callset was joint-called across a
        larger cohort is trimmed first and reported separately.
        """
        parts = []
        if self.n_multiallelic:
            what = ("collapsed (all alleles counted)" if self.multiallelic == "collapse"
                    else "skipped (--multiallelic skip)")
            parts.append(f"{self.n_multiallelic:,} multiallelic record(s) {what}")
        if self.n_alt_trimmed:
            parts.append(f"{self.n_alt_trimmed:,} record(s) lost ALT allele(s) with no reads "
                         f"in this cohort")
        if self.n_monomorphic:
            parts.append(f"{self.n_monomorphic:,} record(s) dropped with no ALT reads at all")
        return "; ".join(parts)


def _usable_alt_indices(alts):
    """Indices of ALTs whose AD column counts as an allele.

    Symbolic ALTs (``<NON_REF>``, ``<*>``) are not alleles and their column is dropped.
    The spanning-deletion placeholder ``*`` is kept: reads carrying a deletion over the
    site are a distinct haplotype, so a sample split between them and a base is mixed
    (moimix counts them too, which is what keeps the two in agreement).
    """
    return [i for i, a in enumerate(alts) if a and not a.startswith("<")]


def _select_alleles(ad, ref, alts, *, trim, snps_only):
    """Choose the AD columns to score at one record.

    ``ad`` is ``[n_samples, 1 + len(alts)]`` with missing depths already zero. Returns
    ``(cols, trimmed)``: the column indices to keep (REF first), or ``None`` when the
    record is dropped; ``trimmed`` says whether an ALT was removed for having no reads.

    With ``trim`` an ALT with no reads in any sample is dropped before anything else is
    decided -- ``bcftools view --trim-alt-alleles`` done on AD rather than GT. A callset
    joint-called across a larger cohort carries ALTs no sample here supports; they add
    nothing to the estimate (zero depth is zero frequency), but they would make a SNP
    site look like an indel site to ``snps_only`` and a biallelic one look multiallelic.
    ``snps_only`` then needs REF and every remaining ALT to be a single base.
    """
    keep = _usable_alt_indices(alts)
    trimmed = False
    if trim:
        present = ad.sum(axis=0) > 0
        kept = [i for i in keep if present[i + 1]]
        trimmed = len(kept) < len(keep)
        keep = kept
    if not keep:
        return None, trimmed
    if snps_only and (len(ref) != 1 or any(len(alts[i]) != 1 for i in keep)):
        return None, trimmed
    return [0] + [i + 1 for i in keep], trimmed


class _Tally:
    """Per-reader bookkeeping shared by the two front-ends."""

    def __init__(self, multiallelic):
        _check_multiallelic(multiallelic)
        self.multiallelic = multiallelic
        self.records, self.sites, self.alleles = [], [], []
        self.n_multi = self.n_alt_trimmed = self.n_monomorphic = 0

    def take(self, ad, cols, trimmed, *, site=None, names=None):
        """Record one site's chosen columns, or account for why it was dropped.

        ``site`` is the site id and ``names`` the allele strings for *all* of ``ad``'s
        columns (REF first); the kept ones are stored in column order.
        """
        if cols is None:
            if trimmed and ad.sum() == ad[:, 0].sum():
                self.n_monomorphic += 1  # every ALT read-less: nothing to score
            return
        if trimmed:
            self.n_alt_trimmed += 1
        if len(cols) > 2:
            self.n_multi += 1
            if self.multiallelic == "skip":
                return
        self.records.append(ad[:, cols])
        self.sites.append(site)
        self.alleles.append(None if names is None else [names[c] for c in cols])

    def finish(self, n_samples):
        return AlleleDepths.from_records(
            self.records, n_samples, n_multiallelic=self.n_multi,
            multiallelic=self.multiallelic, n_alt_trimmed=self.n_alt_trimmed,
            n_monomorphic=self.n_monomorphic, sites=self.sites, alleles=self.alleles)


# --------------------------------------------------------------------------- #
#  Front-ends: AD table  or  VCF/BCF  ->  AlleleDepths                        #
# --------------------------------------------------------------------------- #


def read_ad_table(path, samples, exclude=None, snps_only=False, multiallelic="collapse",
                  trim=True):
    """Read a bcftools-query TSV (``CHROM POS REF ALT`` then one AD per sample).

    The AD column holds the full comma-separated FORMAT/AD (``ref,alt[,alt2...]``), as
    ``bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n'`` writes it. Returns an
    :class:`AlleleDepths` over the sites kept, dropping sites inside any excluded window.

    With ``trim`` (the default) an ALT with no reads in any sample is dropped before the
    site is classified, so a joint-called callset's cohort-wide alleles neither turn a
    SNP site into a non-SNP one nor count as multiallelic here. ``snps_only=True`` then
    requires a single-base REF and single-base ALTs. A multiallelic row (more than one
    ALT with reads) is collapsed with every allele's depth counted, or dropped with
    ``multiallelic="skip"``; either way it is counted in ``n_multiallelic``. A per-sample
    AD of ``.`` or an empty field counts as zero depth.
    """
    exclude = exclude or {}
    tally = _Tally(multiallelic)
    n_mismatch = 0
    seen_counts: set[int] = set()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, pos, r, a = f[0], f[1], f[2], f[3]
            if a == ".":
                continue
            alts = a.split(",")
            if _in_excluded(chrom, int(pos), exclude):
                continue
            ads = f[4:]
            if len(ads) != len(samples):
                n_mismatch += 1
                seen_counts.add(len(ads))
                continue
            n_alleles = 1 + len(alts)
            ad = np.zeros((len(samples), n_alleles), dtype=float)
            ok = True
            for j, cell in enumerate(ads):
                parts = cell.split(",")
                if len(parts) < 2 or parts[0] in (".", ""):
                    continue  # missing -> zero depth
                if len(parts) > n_alleles:
                    ok = False
                    break
                try:
                    ad[j, :len(parts)] = [0 if x in (".", "") else int(x) for x in parts]
                except ValueError:
                    ok = False
                    break
            if ok:
                tally.take(ad, *_select_alleles(ad, r, alts, trim=trim, snps_only=snps_only),
                           site=f"{chrom}:{pos}", names=[r] + alts)
    if not tally.records and n_mismatch:
        raise ValueError(
            f"AD table has {sorted(seen_counts)} value column(s) per site but "
            f"{len(samples)} samples were given — all {n_mismatch} sites dropped. The AD "
            f"columns and sample list are misaligned (e.g. a multi-@RG-SM BAM adds columns; "
            f"run mpileup with --ignore-RG).")
    return tally.finish(len(samples))


def read_ad_vcf(path, exclude=None, snps_only=False, multiallelic="collapse", trim=True):
    """Read per-sample AD from a VCF/BCF.

    Returns ``(samples, depths)``: the sample id list plus an :class:`AlleleDepths` over
    every record carrying ``AD`` and at least one ALT with reads. Pass the **unsplit**
    callset: a multiallelic record is collapsed with every allele's depth counted (the
    default), or dropped with ``multiallelic="skip"``; ``bcftools norm -m-`` output cannot
    be un-split, and each split record has already discarded the other alleles' reads.

    With ``trim`` (the default) an ALT with no reads in any sample is dropped before the
    site is classified -- see :func:`read_ad_table`. ``snps_only=True`` keeps only
    records whose REF and every remaining ALT are single bases. A missing per-sample AD
    counts as zero depth.
    """
    from cyvcf2 import VCF

    exclude = exclude or {}
    tally = _Tally(multiallelic)
    vcf = VCF(path)
    samples = list(vcf.samples)
    try:
        for v in vcf:
            alts = v.ALT
            if not alts:
                continue
            if exclude and _in_excluded(v.CHROM, v.POS, exclude):
                continue
            ad = v.format("AD")
            if ad is None or ad.shape[1] < 1 + len(alts):
                continue
            ad = np.where(ad < 0, 0, ad).astype(float)  # cyvcf2 missing sentinel -> 0 depth
            tally.take(ad, *_select_alleles(ad, v.REF, alts, trim=trim, snps_only=snps_only),
                       site=f"{v.CHROM}:{v.POS}", names=[v.REF] + list(alts))
    finally:
        vcf.close()
    return samples, tally.finish(len(samples))


def read_allele_table(path, *, sample_col="library_sample_name", locus_col="target_name",
                      allele_col="seq", reads_col="reads"):
    """Read a long-format allele table (one row per sample, locus and allele) into depths.

    This is the shape amplicon pipelines write -- the defaults are MAD4HATTER's columns --
    and it turns each locus into one record whose alleles are the haplotypes seen anywhere
    in the cohort and whose per-sample depths are the read counts, so microhaplotypes go
    through the same estimator as SNPs. A sample with no row at a locus has zero depth
    there. Rows with the same sample, locus and allele are summed. Returns
    ``(samples, depths)`` with samples in sorted order.

    There is no REF here: allele 0 of every record is whichever haplotype came first, so
    ``compute_fws(min_alt_samples=...)`` is not meaningful on these depths. Trimming and
    the SNP rule do not apply either -- every allele present has reads by construction.
    ``n_multiallelic`` counts loci with more than one allele.
    """
    import pandas as pd

    df = pd.read_csv(path, sep="\t", usecols=[sample_col, locus_col, allele_col, reads_col],
                     dtype={sample_col: str, locus_col: str, allele_col: str})
    df = df.dropna(subset=[sample_col, locus_col, allele_col])
    df[reads_col] = pd.to_numeric(df[reads_col], errors="coerce").fillna(0)
    samples = sorted(df[sample_col].unique())
    sidx = {s: i for i, s in enumerate(samples)}
    records, sites, names, n_multi = [], [], [], 0
    for locus, g in df.groupby(locus_col, sort=True):
        alleles = {a: k for k, a in enumerate(g[allele_col].unique())}
        sites.append(str(locus))
        names.append(list(alleles))
        ad = np.zeros((len(samples), len(alleles)), dtype=float)
        np.add.at(ad, (g[sample_col].map(sidx).values, g[allele_col].map(alleles).values),
                  g[reads_col].values)
        if len(alleles) > 1:
            n_multi += 1
        records.append(ad)
    depths = AlleleDepths.from_records(records, len(samples), n_multiallelic=n_multi,
                                       multiallelic="collapse", sites=sites, alleles=names)
    return samples, depths


def read_pop_freqs(path, *, locus_col="locus", allele_col="allele", freq_col="freq"):
    """Read population allele frequencies from a long-format TSV into ``{site: {allele: freq}}``.

    One row per locus and allele. For VCF / AD-table input the locus is ``CHROM:POS`` and
    the allele its REF or ALT string; for an allele table it is the locus name and the
    haplotype. This is the format :func:`write_pop_freqs` writes, so frequencies computed
    on a reference cohort can be reused. Frequencies are renormalised per locus to sum to
    1, so counts work too.
    """
    import pandas as pd

    df = pd.read_csv(path, sep="\t", usecols=[locus_col, allele_col, freq_col],
                     dtype={locus_col: str, allele_col: str})
    df[freq_col] = pd.to_numeric(df[freq_col], errors="coerce").fillna(0.0)
    out: dict[str, dict[str, float]] = {}
    for locus, g in df.groupby(locus_col, sort=False):
        tot = float(g[freq_col].sum())
        if tot <= 0:
            continue
        out[str(locus)] = {a: float(f) / tot for a, f in zip(g[allele_col], g[freq_col])}
    return out


def write_pop_freqs(freqs, path, *, locus_col="locus", allele_col="allele", freq_col="freq"):
    """Write ``{site: {allele: freq}}`` (see :meth:`AlleleDepths.population_freqs`) as a TSV."""
    with open(path, "w") as fh:
        fh.write(f"{locus_col}\t{allele_col}\t{freq_col}\n")
        for site, d in freqs.items():
            for a, f in d.items():
                fh.write(f"{site}\t{a}\t{f:.6g}\n")


# --------------------------------------------------------------------------- #
#  The estimator                                                              #
# --------------------------------------------------------------------------- #


def compute_fws(depths, alt=None, *, estimator="regression", min_depth=0, n_bins=10,
                min_alt_samples=0, pop_freqs=None):
    """Compute per-sample Fws from an :class:`AlleleDepths`.

    ``compute_fws(ref, alt, ...)`` with two ``[n_sites, n_samples]`` biallelic depth
    matrices is accepted too and is exactly the two-allele case.

    Returns ``(fws, n_sites)`` arrays of length ``n_samples`` (``fws`` is NaN for a
    sample with no usable sites; ``n_sites`` is how many sites it contributed).

    ``estimator`` selects ``"regression"`` (matches ``moimix::getFws``; the default) or
    ``"ratio"``. ``min_depth`` drops per-sample sites below that read depth;
    ``min_alt_samples`` keeps only sites where a non-reference allele is seen in at least
    that many samples. moimix parity uses ``estimator="regression", min_depth=0,
    min_alt_samples=0``.

    ``pop_freqs`` (``{site: {allele: freq}}``, see :func:`read_pop_freqs`) replaces the
    cohort's own pooled read fractions as the population allele frequencies -- for
    scoring a few samples against a reference population, or the same population at a
    different time. Only the population side (``Hs`` and the binning variable) changes;
    the within-sample side never needs frequencies. A site with no entry is dropped, and
    ``compute_fws.last_pop_freq_misses`` says how many were. The population is exactly
    the alleles listed for a site, so an allele seen here but absent there counts as
    frequency 0 and one listed there but unseen here still contributes to ``Hs``.

    ``n_bins=0`` skips the MAF binning and works per site: the regression becomes the
    through-origin slope of every usable site's ``Hw`` on its ``Hs``
    (``Fws = 1 - Σ Hs·Hw / Σ Hs²``), the ratio ``1 - Σ Hw / Σ Hs``. Use it for
    microhaplotype loci, where the major allele is often below 50% and the bins over
    [0, 0.5] stop describing the sites. It is a different estimator from the binned one
    (every site is weighted by its own Hs², rather than each bin by its mean), so on a SNP
    callset the two can differ by a tenth; do not carry a threshold from one to the other.
    """
    if alt is not None:
        depths = AlleleDepths.from_ref_alt(depths, alt)
    if not isinstance(depths, AlleleDepths):
        raise TypeError("compute_fws takes an AlleleDepths, or (ref, alt) matrices")
    if depths.n_sites == 0:
        return np.full(depths.n_samples, np.nan), np.zeros(depths.n_samples, dtype=int)
    depth = depths.depth
    n_samples = depths.n_samples

    if pop_freqs is not None:
        # population from outside: Hs and the minor-allele fraction per site, from the
        # supplied frequencies; a site without an entry is NaN and so drops out below
        if depths.sites is None:
            raise ValueError("pop_freqs needs depths with site ids (a reader built them)")
        Hs = np.full(depths.n_sites, np.nan)
        maf = np.full(depths.n_sites, np.nan)
        for i, site in enumerate(depths.sites):
            f = pop_freqs.get(site)
            if not f:
                continue
            v = np.asarray(list(f.values()), dtype=float)
            v = v / v.sum()
            Hs[i] = 1.0 - float((v * v).sum())
            maf[i] = 1.0 - float(v.max())
        compute_fws.last_pop_freq_misses = int(np.isnan(maf).sum())
    else:
        # population: allele frequencies over the cohort's pooled reads
        tot_dp = depths.pop.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            freq = np.where(tot_dp[:, None] > 0, depths.pop / tot_dp[:, None], np.nan)
            # minor-allele fraction as a ratio of counts (moimix: min(coverage / sum(coverage)))
            maf = np.where(tot_dp > 0, (tot_dp - depths.pop.max(axis=1)) / tot_dp, np.nan)
        Hs = 1.0 - (freq * freq).sum(axis=1)
        compute_fws.last_pop_freq_misses = 0
    alt_present = (depths.nonref > 0).sum(axis=1)
    # within-sample: Hw = 1 - Σ_k q_k² = 1 - Σ_k ad_k² / depth²
    with np.errstate(invalid="ignore", divide="ignore"):
        Hw = np.where(depth > 0, 1.0 - depths.sumsq / (depth * depth), np.nan)

    fws = np.full(n_samples, np.nan)
    n_info = np.zeros(n_samples, dtype=int)

    if n_bins == 0:
        # per-site, no binning: the same through-origin regression (or ratio) with every
        # usable site as its own point
        site_ok = np.isfinite(maf) & np.isfinite(Hs) & (alt_present >= min_alt_samples)
        if estimator == "ratio":
            site_ok &= Hs > 0
        elif estimator != "regression":
            raise ValueError(f"unknown estimator {estimator!r} (use 'regression' or 'ratio')")
        usable = site_ok[:, None] & (depth > 0) & (depth >= min_depth) & np.isfinite(Hw)
        x = np.where(usable, Hs[:, None], 0.0)
        y = np.where(usable, Hw, 0.0)
        y = np.where(np.isfinite(y), y, 0.0)
        denom = (x * x).sum(axis=0) if estimator == "regression" else x.sum(axis=0)
        numer = (x * y).sum(axis=0) if estimator == "regression" else y.sum(axis=0)
        ok = denom > 0
        fws[ok] = 1.0 - numer[ok] / denom[ok]
        n_info[:] = usable.sum(axis=0)
        n_info[~ok] = 0
        return fws, n_info

    # Bins of width 0.5 / n_bins over [0, 0.5] -- moimix's grid, built for a biallelic
    # site whose minor-allele fraction cannot exceed 0.5. A site with k alleles has
    # `maf = 1 - max(p)` up to 1 - 1/k, and used to fall off the top of that grid: the
    # regression estimator gave every such site one shared overflow bin (findInterval's
    # n_bins + 1), the ratio estimator clipped them into the 0.45-0.5 bin beside biallelic
    # sites of half their heterozygosity. Neither is a bin of *similar* sites, which is what
    # the binning is for. So keep the biallelic grid exactly as it was -- nothing biallelic
    # moves -- and extend it upward in the same width as far as the data reach.
    width = 0.5 / n_bins
    maf_top = np.nanmax(maf) if np.isfinite(maf).any() else 0.5
    n_extra = int(np.ceil(max(0.0, maf_top - 0.5) / width - 1e-9)) if maf_top > 0.5 else 0
    edges = np.linspace(0, 0.5 + n_extra * width, n_bins + n_extra + 1)
    n_bins_eff = n_bins + n_extra

    if estimator == "regression":
        # moimix::getFws — 10 MAF bins via findInterval, global per-bin population-het
        # means, Fws = 1 - slope of a through-origin regression of the per-sample binned
        # het means on those population means.
        site_ok = np.isfinite(maf) & (alt_present >= min_alt_samples)
        if not site_ok.any():
            return fws, n_info
        bin_idx = np.searchsorted(edges, maf, side="right")  # findInterval: 1..n_bins+1
        bins = np.unique(bin_idx[site_ok])
        xbar = {}
        for b in bins:
            sel = site_ok & (bin_idx == b)
            xbar[b] = np.nanmean(Hs[sel]) if sel.any() else np.nan
        for s in range(n_samples):
            usable = site_ok & (depth[:, s] > 0) & (depth[:, s] >= min_depth)
            if not usable.any():
                continue
            xs, ys = [], []
            for b in bins:
                sel = usable & (bin_idx == b)
                if not sel.any():
                    continue
                y = np.nanmean(Hw[sel, s])
                x = xbar[b]
                if np.isfinite(x) and np.isfinite(y):
                    xs.append(x)
                    ys.append(y)
            xs = np.asarray(xs)
            ys = np.asarray(ys)
            denom = float((xs * xs).sum())
            if denom > 0:
                fws[s] = 1.0 - float((xs * ys).sum()) / denom
                n_info[s] = int(usable.sum())
        return fws, n_info

    if estimator == "ratio":
        # Fws = 1 - Σ_bins mean(Hw) / Σ_bins mean(Hs), over polymorphic sites, with
        # per-sample bin means.
        site_ok = (Hs > 0) & (alt_present >= min_alt_samples) & np.isfinite(maf)
        bin_idx = np.clip(np.digitize(maf, edges[1:-1]), 0, n_bins_eff - 1)
        for s in range(n_samples):
            usable = site_ok & (depth[:, s] >= min_depth)
            if not usable.any():
                continue
            sum_hw = sum_hs = 0.0
            for b in range(n_bins_eff):
                sel = usable & (bin_idx == b)
                if sel.any():
                    sum_hw += np.nanmean(Hw[sel, s])
                    sum_hs += np.nanmean(Hs[sel])
            if sum_hs > 0:
                fws[s] = 1.0 - sum_hw / sum_hs
                n_info[s] = int(usable.sum())
        return fws, n_info

    raise ValueError(f"unknown estimator {estimator!r} (use 'regression' or 'ratio')")


def fws_table(path, *, fws_min=0.95, estimator="regression", min_depth=0, n_bins=10,
              min_alt_samples=0, snps_only=True, multiallelic="collapse", trim=True,
              exclude_call_regions=None):
    """Score every sample in a callset and say which ones a ``fws_min`` cut would keep.

    Returns ``(rows, n_sites)``: one row per sample with ``sample``, ``fws`` (``None`` where
    it could not be computed), ``n_sites``, ``monoclonal`` and ``dropped``, in file order.

    A sample with no usable sites cannot be scored, and an unscored sample is **not** a
    monoclonal one -- it is a sample nothing is known about, so it is dropped and counted
    apart from the polyclonal ones. Silently keeping it would put exactly the samples this
    step exists to exclude back in the output.

    ``multiallelic`` is ``"collapse"`` (count every allele's depth at a multiallelic site,
    the default) or ``"skip"`` (drop such sites); ``trim`` drops ALTs no sample here has
    reads for before a site is classified -- see :func:`read_ad_vcf`.
    """
    exclude = load_exclude_regions(exclude_call_regions) if exclude_call_regions else {}
    samples, depths = read_ad_vcf(path, exclude, snps_only=snps_only,
                                  multiallelic=multiallelic, trim=trim)
    if depths.n_sites == 0:
        raise SystemExit("fws_filter: no usable " + ("SNP" if snps_only else "variant")
                         + " sites with AD in the input")
    if depths.multiallelic_note():
        detail(f"       {depths.multiallelic_note()}")
    fws, n_info = compute_fws(depths, estimator=estimator, min_depth=min_depth,
                              n_bins=n_bins, min_alt_samples=min_alt_samples)
    rows = []
    for s, f, n in zip(samples, fws, n_info):
        scored = bool(np.isfinite(f))
        mono = scored and float(f) >= fws_min
        rows.append({"sample": s, "fws": float(f) if scored else None,
                     "n_sites": int(n), "monoclonal": mono, "dropped": not mono})
    return rows, depths.n_sites


def write_fws_table(rows, path):
    """Write what :func:`fws_table` decided, one row per sample."""
    import csv as _csv

    with open(path, "w", newline="") as fh:
        w = _csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "fws", "n_sites", "monoclonal", "dropped"])
        for r in rows:
            w.writerow([r["sample"], "" if r["fws"] is None else f"{r['fws']:.6f}",
                        r["n_sites"], r["monoclonal"], r["dropped"]])


def fws_filter(inp, out, *, fws_min=0.95, estimator="regression", min_depth=0, n_bins=10,
               min_alt_samples=0, snps_only=True, multiallelic="collapse", trim=True,
               exclude_call_regions=None, fws_table_path=None, dropped_samples_path=None):
    """Keep only samples with Fws >= ``fws_min`` -- the monoclonal infections. Drops samples.

    This is an analysis choice rather than a QC rule, and it is the one step in the chain
    that changes *which infections* the callset describes, so it is off in the default
    config and has to be asked for.

    **It removes no variants.** Sites the remaining samples no longer support are left in
    place: dropping samples changes every allele frequency, so a site that cleared a MAF or
    missingness bar with the whole cohort may not clear it with this one. That is a real
    consequence and not one to bury inside a sample filter -- put ``maf_filter`` and
    ``locus_missingness_filter`` after this step to re-apply them to the survivors.
    ``AC``/``AN``/``AF`` are refreshed here, so those steps read the new frequencies.

    Fws is measured against the cohort's own allele frequencies, so it wants the callset
    the rest of the chain has already cleaned -- run it at the end, not as an entry gate.
    Re-genotyping upstream is what removes the minor-allele noise Fws would otherwise read
    as within-host diversity.

    ``fws_table_path`` writes the per-sample scores the decision was made from; a sample
    that vanished from a cohort is otherwise just a name in a log.

    Returns the list of dropped sample names.
    """
    from .bcftools import out_flag, q, require, sh

    require("bcftools")
    rows, n_sites = fws_table(inp, fws_min=fws_min, estimator=estimator,
                              min_depth=min_depth, n_bins=n_bins,
                              min_alt_samples=min_alt_samples, snps_only=snps_only,
                              multiallelic=multiallelic, trim=trim,
                              exclude_call_regions=exclude_call_regions)
    dropped = sorted(r["sample"] for r in rows if r["dropped"])
    unscored = sorted(r["sample"] for r in rows if r["fws"] is None)
    if len(dropped) == len(rows):
        raise SystemExit(
            f"fws_filter: Fws >= {fws_min:g} keeps no samples of {len(rows)} "
            f"(scored over {n_sites:,} site(s)). Lower the threshold, or check that this "
            f"callset is filtered and re-genotyped -- unfiltered calls read as within-host "
            f"diversity and push every sample's Fws down.")

    if fws_table_path:
        write_fws_table(rows, fws_table_path)

    # the borderline ones said out loud, the way sample_coverage_filter does: a drop that
    # missed by a hair is the one worth seeing without opening the table
    for r in sorted(rows, key=lambda r: (r["fws"] is not None, r["fws"] or 0)):
        near = r["fws"] is not None and abs(r["fws"] - fws_min) <= 0.05
        if r["dropped"] or near:
            score = "unscored" if r["fws"] is None else f"{r['fws']:.4f}"
            margin = "" if r["fws"] is None else f" ({r['fws'] - fws_min:+.4f})"
            detail(f"       {r['sample']}\tFws {score}\t{r['n_sites']:,} sites"
                  f"\t{'DROPPED' if r['dropped'] else 'kept'}{margin}")

    if dropped_samples_path:
        with open(dropped_samples_path, "w") as fh:
            fh.write("\n".join(dropped) + ("\n" if dropped else ""))

    fmt = out_flag(out)
    if dropped:
        from .vcf_filters import _write_tmp_list
        drop_arg = dropped_samples_path or _write_tmp_list(dropped)
        view = f"bcftools view -S ^{q(drop_arg)} {q(inp)} -Ou"
    else:
        view = f"bcftools view {q(inp)} -Ou"
    sh(f"{view} | bcftools +fill-tags -O{fmt} -o {q(out)} -- -t AC,AN,AF",
       tools=("bcftools",))
    if unscored:
        say(f"     {len(unscored)} sample(s) could not be scored and were dropped: "
              + ", ".join(unscored))
    return dropped
