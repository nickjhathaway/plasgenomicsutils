"""Fws within-host diversity statistic (Manske 2012; the moimix::getFws estimator).

Fws measures how much of a sample's genetic diversity is captured within the host
versus the population — a monoclonal infection has Fws ~ 1, a polyclonal one less.
It needs a cohort (to estimate population allele frequencies) and per-sample allele
depths (AD) at biallelic SNP sites:

  * population alt frequency at site i:  p_i = Σ_s alt_si / Σ_s (ref_si + alt_si)
  * population heterozygosity:           Hs_i = 2 p_i (1 - p_i)
  * within-sample alt frequency:         q_si = alt_si / (ref_si + alt_si)
  * within-sample heterozygosity:        Hw_si = 2 q_si (1 - q_si)
  * sites are binned by population MAF into ``n_bins`` equal bins over [0, 0.5];
    per bin the mean Hw (per sample) and mean Hs are taken.

Two estimators combine those binned means into Fws — see :func:`compute_fws`:

  * ``"regression"`` (default): Fws = 1 - β, where β is the slope of a regression of
    the binned sample-het means on the binned population-het means, forced through
    the origin. This matches ``moimix::getFws``.
  * ``"ratio"``: Fws = 1 - Σ_bins mean(Hw) / Σ_bins mean(Hs) — a simpler ratio of the
    summed binned means.

The two agree in spirit but not to the digit (the regression weights bins by the
squared population het), so pick deliberately and don't mix a threshold tuned on one
with values from the other.

AD is read from either a bcftools-query **AD table** (:func:`read_ad_table`) or a
**VCF/BCF** (:func:`read_ad_vcf`), both feeding the same :func:`compute_fws`.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .reporting import detail, say

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
#  Front-ends: AD table  or  VCF/BCF  ->  (ref, alt) depth matrices           #
# --------------------------------------------------------------------------- #


def read_ad_table(path, samples, exclude=None, snps_only=False):
    """Read a bcftools-query TSV (``CHROM POS REF ALT`` then one ``ref,alt`` AD per sample).

    Returns ``(ref, alt)`` float arrays of shape ``[n_sites, n_samples]`` over biallelic
    sites, dropping sites inside any excluded window. Every row with a real ALT is kept
    (matching :func:`read_ad_vcf` and moimix, which use AD regardless of allele string);
    ``snps_only=True`` additionally requires single-base REF/ALT. A per-sample AD of ``.``
    or an empty field counts as zero depth.
    """
    exclude = exclude or {}
    ref_rows, alt_rows = [], []
    n_mismatch = 0
    seen_counts: set[int] = set()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, pos, r, a = f[0], f[1], f[2], f[3]
            if a in (".", "<*>"):
                continue
            if snps_only and (len(r) != 1 or len(a) != 1):
                continue
            if _in_excluded(chrom, int(pos), exclude):
                continue
            ads = f[4:]
            if len(ads) != len(samples):
                n_mismatch += 1
                seen_counts.add(len(ads))
                continue
            refs, alts, ok = [], [], True
            for ad in ads:
                parts = ad.split(",")
                if len(parts) < 2 or parts[0] in (".", ""):
                    refs.append(0)
                    alts.append(0)
                else:
                    try:
                        refs.append(int(parts[0]))
                        alts.append(int(parts[1]))
                    except ValueError:
                        ok = False
                        break
            if ok:
                ref_rows.append(refs)
                alt_rows.append(alts)
    if not ref_rows and n_mismatch:
        raise ValueError(
            f"AD table has {sorted(seen_counts)} value column(s) per site but "
            f"{len(samples)} samples were given — all {n_mismatch} sites dropped. The AD "
            f"columns and sample list are misaligned (e.g. a multi-@RG-SM BAM adds columns; "
            f"run mpileup with --ignore-RG).")
    return np.array(ref_rows, dtype=float), np.array(alt_rows, dtype=float)


def read_ad_vcf(path, exclude=None, snps_only=False):
    """Read per-sample AD from a VCF/BCF over biallelic sites.

    Returns ``(samples, ref, alt)``: the sample id list plus float ``[n_sites, n_samples]``
    ref/alt depth arrays, taken from the first two AD columns (ref, alt) of every biallelic
    record — matching ``moimix::getFws``, which uses AD regardless of allele string. Pass a
    biallelic callset (``bcftools norm -m-``); this reads whatever survives, so filter to the
    site set you want (e.g. SNPs) upstream. ``snps_only=True`` additionally keeps only
    single-base REF/ALT records. A missing per-sample AD counts as zero depth.
    """
    from cyvcf2 import VCF

    exclude = exclude or {}
    vcf = VCF(path)
    samples = list(vcf.samples)
    ref_rows, alt_rows = [], []
    try:
        for v in vcf:
            if len(v.ALT) != 1:
                continue  # biallelic only (split multiallelics first)
            if snps_only and (len(v.REF) != 1 or len(v.ALT[0]) != 1):
                continue
            if exclude and _in_excluded(v.CHROM, v.POS, exclude):
                continue
            ad = v.format("AD")
            if ad is None or ad.shape[1] < 2:
                continue
            ad = np.where(ad < 0, 0, ad).astype(float)  # cyvcf2 missing sentinel -> 0 depth
            ref_rows.append(ad[:, 0])
            alt_rows.append(ad[:, 1])
    finally:
        vcf.close()
    ref = np.array(ref_rows, dtype=float) if ref_rows else np.empty((0, len(samples)))
    alt = np.array(alt_rows, dtype=float) if alt_rows else np.empty((0, len(samples)))
    return samples, ref, alt


# --------------------------------------------------------------------------- #
#  The estimator                                                              #
# --------------------------------------------------------------------------- #


def compute_fws(ref, alt, *, estimator="regression", min_depth=0, n_bins=10,
                min_alt_samples=0):
    """Compute per-sample Fws from ``[n_sites, n_samples]`` ref/alt depth matrices.

    Returns ``(fws, n_sites)`` arrays of length ``n_samples`` (``fws`` is NaN for a
    sample with no usable sites; ``n_sites`` is how many sites it contributed).

    ``estimator`` selects ``"regression"`` (matches ``moimix::getFws``; the default) or
    ``"ratio"``. ``min_depth`` drops per-sample sites below that read depth;
    ``min_alt_samples`` keeps only sites where the alt is seen in at least that many
    samples. moimix parity uses ``estimator="regression", min_depth=0,
    min_alt_samples=0``.
    """
    ref = np.asarray(ref, dtype=float)
    alt = np.asarray(alt, dtype=float)
    if ref.size == 0:
        return np.full(0, np.nan), np.zeros(0, dtype=int)
    depth = ref + alt
    _n_sites, n_samples = ref.shape

    tot_alt = alt.sum(axis=1)
    tot_dp = depth.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(tot_dp > 0, tot_alt / tot_dp, np.nan)
    maf = np.minimum(p, 1 - p)
    Hs = 2 * p * (1 - p)
    alt_present = (alt > 0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(depth > 0, alt / depth, np.nan)
    Hw = 2 * q * (1 - q)

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
              min_alt_samples=0, snps_only=True, exclude_call_regions=None):
    """Score every sample in a callset and say which ones a ``fws_min`` cut would keep.

    Returns ``(rows, n_sites)``: one row per sample with ``sample``, ``fws`` (``None`` where
    it could not be computed), ``n_sites``, ``monoclonal`` and ``dropped``, in file order.

    A sample with no usable sites cannot be scored, and an unscored sample is **not** a
    monoclonal one -- it is a sample nothing is known about, so it is dropped and counted
    apart from the polyclonal ones. Silently keeping it would put exactly the samples this
    step exists to exclude back in the output.
    """
    import numpy as np

    exclude = load_exclude_regions(exclude_call_regions) if exclude_call_regions else {}
    samples, ref, alt = read_ad_vcf(path, exclude, snps_only=snps_only)
    if ref.size == 0:
        raise SystemExit("fws_filter: no usable biallelic SNP sites in the input")
    fws, n_info = compute_fws(ref, alt, estimator=estimator, min_depth=min_depth,
                              n_bins=n_bins, min_alt_samples=min_alt_samples)
    rows = []
    for s, f, n in zip(samples, fws, n_info):
        scored = bool(np.isfinite(f))
        mono = scored and float(f) >= fws_min
        rows.append({"sample": s, "fws": float(f) if scored else None,
                     "n_sites": int(n), "monoclonal": mono, "dropped": not mono})
    return rows, int(ref.shape[0])


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
               min_alt_samples=0, snps_only=True, exclude_call_regions=None,
               fws_table_path=None, dropped_samples_path=None):
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
