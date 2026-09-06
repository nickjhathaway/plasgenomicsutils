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

AD is read from either a bcftools-query **AD table** (:func:`read_ad_table`) or a
**VCF/BCF** (:func:`read_ad_vcf`) into an :class:`AlleleDepths`, which feeds
:func:`compute_fws`.
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
        ``multiallelic``, ``n_alt_trimmed``, ``n_monomorphic``).
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
        self.records = []
        self.n_multi = self.n_alt_trimmed = self.n_monomorphic = 0

    def take(self, ad, cols, trimmed):
        """Record one site's chosen columns, or account for why it was dropped."""
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

    def finish(self, n_samples):
        return AlleleDepths.from_records(
            self.records, n_samples, n_multiallelic=self.n_multi,
            multiallelic=self.multiallelic, n_alt_trimmed=self.n_alt_trimmed,
            n_monomorphic=self.n_monomorphic)


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
                tally.take(ad, *_select_alleles(ad, r, alts, trim=trim, snps_only=snps_only))
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
            tally.take(ad, *_select_alleles(ad, v.REF, alts, trim=trim, snps_only=snps_only))
    finally:
        vcf.close()
    return samples, tally.finish(len(samples))


# --------------------------------------------------------------------------- #
#  The estimator                                                              #
# --------------------------------------------------------------------------- #


def compute_fws(depths, alt=None, *, estimator="regression", min_depth=0, n_bins=10,
                min_alt_samples=0):
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
    """
    if alt is not None:
        depths = AlleleDepths.from_ref_alt(depths, alt)
    if not isinstance(depths, AlleleDepths):
        raise TypeError("compute_fws takes an AlleleDepths, or (ref, alt) matrices")
    if depths.n_sites == 0:
        return np.full(depths.n_samples, np.nan), np.zeros(depths.n_samples, dtype=int)
    depth = depths.depth
    n_samples = depths.n_samples

    # population: allele frequencies over the cohort's pooled reads
    tot_dp = depths.pop.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        freq = np.where(tot_dp[:, None] > 0, depths.pop / tot_dp[:, None], np.nan)
        # minor-allele fraction as a ratio of counts (moimix: min(coverage / sum(coverage)))
        maf = np.where(tot_dp > 0, (tot_dp - depths.pop.max(axis=1)) / tot_dp, np.nan)
    Hs = 1.0 - (freq * freq).sum(axis=1)
    alt_present = (depths.nonref > 0).sum(axis=1)
    # within-sample: Hw = 1 - Σ_k q_k² = 1 - Σ_k ad_k² / depth²
    with np.errstate(invalid="ignore", divide="ignore"):
        Hw = np.where(depth > 0, 1.0 - depths.sumsq / (depth * depth), np.nan)

    edges = np.linspace(0, 0.5, n_bins + 1)
    fws = np.full(n_samples, np.nan)
    n_info = np.zeros(n_samples, dtype=int)

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
        bin_idx = np.clip(np.digitize(maf, edges[1:-1]), 0, n_bins - 1)
        for s in range(n_samples):
            usable = site_ok & (depth[:, s] >= min_depth)
            if not usable.any():
                continue
            sum_hw = sum_hs = 0.0
            for b in range(n_bins):
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
