"""What the callset looks like at the end: two report tables that change nothing.

``sample_summary`` is the per-sample view -- the coverage columns ``sample_coverage_filter``
decides on and the Fws columns ``fws_filter`` decides on, side by side, with no decision
made. ``variant_summary`` is the per-record view: every record classed, and within each
class broken down by how many ALT alleles it carries, as counts and as fractions.

Both exist because a chain that ends in ``maf_filter`` tells you how many variants
survived and nothing about the callset those variants make -- which samples are thin, which
are mixed, how much of it is biallelic. The filters that *would* say so are the ones most
often switched off (``sample_coverage_filter``, ``fws_filter``), and they say it only for
the samples they drop.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile

from .bcftools import classify_record, q, require, sh
from .reporting import say

#: How a record with this many ALT alleles is named in the summary.
ALLELE_LABELS = {1: "biallelic", 2: "triallelic", 3: "tetra-allelic"}


def _allele_label(n_alt: int) -> str:
    return ALLELE_LABELS.get(n_alt, f"{n_alt + 1}-allelic")


def variant_summary_table(path: str) -> list[dict]:
    """Records by class and by ALT-allele count, with their share of the class and of all.

    One row per ``(class, n_alt)`` present, classes in :data:`VARIANT_TYPES` order and
    allele counts ascending; a ``(class, "all")`` row heads each class. ``frac_total`` is the
    share of every record in the file, ``frac_class`` the share within the class, so
    "100 biallelic SNPs, 0.50 of the callset, 0.98 of the SNPs" is one row.
    """
    from .bcftools import VARIANT_TYPES

    require("bcftools")
    proc = subprocess.run(f"bcftools query -f '%REF\\t%ALT\\n' {q(path)}", shell=True,
                          executable="/bin/bash", stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: could not read variants from {path}: {proc.stderr.strip()}")
    counts: dict[tuple[str, int], int] = {}
    total = 0
    for line in proc.stdout.splitlines():
        ref, _, alt = line.partition("\t")
        cls = classify_record(ref, alt)
        n_alt = 0 if alt in (".", "") else alt.count(",") + 1
        counts[(cls, n_alt)] = counts.get((cls, n_alt), 0) + 1
        total += 1
    rows = []
    for cls in VARIANT_TYPES:
        per = {k[1]: v for k, v in counts.items() if k[0] == cls}
        if not per:
            continue
        n_cls = sum(per.values())
        rows.append({"class": cls, "n_alt": "all", "alleles": "all", "count": n_cls,
                     "frac_total": round(n_cls / total, 4), "frac_class": 1.0})
        for n_alt in sorted(per):
            rows.append({"class": cls, "n_alt": n_alt, "alleles": _allele_label(n_alt),
                         "count": per[n_alt], "frac_total": round(per[n_alt] / total, 4),
                         "frac_class": round(per[n_alt] / n_cls, 4)})
    rows.append({"class": "total", "n_alt": "all", "alleles": "all", "count": total,
                 "frac_total": 1.0 if total else 0.0, "frac_class": 1.0 if total else 0.0})
    return rows


VARIANT_SUMMARY_COLUMNS = ["class", "n_alt", "alleles", "count", "frac_total", "frac_class"]

#: The frequency marks the spectrum is reported at. These are the floors a downstream
#: analysis is likely to apply, so the table answers "how much of this callset would
#: survive if I filtered at X" without re-reading the callset.
MAF_MARKS = (0.01, 0.02, 0.05, 0.10)

MAF_SPECTRUM_COLUMNS = ["group", "n_samples", "n_records", "maf_min", "at_or_above",
                        "frac_at_or_above", "in_band_below_next"]


def maf_spectrum_table(path: str, *, marks=MAF_MARKS, meta: str | None = None,
                       group_col: str | None = None, sample_col: str = "sample") -> list[dict]:
    """How many records sit at or above each frequency mark, overall and per group.

    Run on the callset **before** the frequency filter, this says where the bulk of the
    variation actually is -- and therefore what a 1%, 2%, 5% or 10% floor would cost. The
    site frequency spectrum of a *P. falciparum* cohort is steep enough that the answer is
    rarely the one a round number suggests: on a 249-sample callset, 39% of records clear
    2% and 62% clear 1%, so the choice of floor moves more of the panel than the two
    thresholds look like they should.

    ``in_band_below_next`` is the count between this mark and the next one up, which is what
    a change of floor actually gains or loses.

    With ``meta`` + ``group_col`` the spectrum is computed **per group as well**, because a
    grouped frequency floor (``maf_filter --meta --group-col``) keeps a record when any one
    group clears the bar: an allele at 3% in one country and absent in another is 1.5%
    pooled, and which number matters depends on how the filter is being run.
    """
    require("bcftools")
    groups: dict[str, list[str]] = {}
    if meta and group_col:
        import csv as _csv

        from ..utils.small_utils import Utils
        from .vcf_filters import _vcf_samples

        present = _vcf_samples(path)
        with open(meta) as fh:
            reader = _csv.DictReader(fh, delimiter="\t")
            fields = reader.fieldnames or []
            s_col = Utils.resolve_column(fields, sample_col, source=f"metadata ({meta})")
            g_col = Utils.resolve_column(fields, group_col, source=f"metadata ({meta})")
            for row in reader:
                if row[s_col] in present and row[g_col]:
                    groups.setdefault(row[g_col], []).append(row[s_col])
    rows = []
    rows += _spectrum_rows(path, "ALL", None, marks)
    for g, names in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        rows += _spectrum_rows(path, g, names, marks)
    return rows


def _spectrum_rows(path: str, label: str, samples: list[str] | None, marks) -> list[dict]:
    """The spectrum for one group (or the whole callset when ``samples`` is None)."""
    tmp = None
    try:
        src = path
        if samples:
            tmp = tempfile.mkdtemp(prefix="maf_spectrum_")
            sfile = os.path.join(tmp, "s.txt")
            with open(sfile, "w") as fh:
                fh.write("\n".join(samples) + "\n")
            src = os.path.join(tmp, "g.bcf")
            sh(f"bcftools view -S {q(sfile)} --force-samples {q(path)} -Ob -o {q(src)}",
               tools=("bcftools",))
        proc = subprocess.run(
            f"bcftools +fill-tags {q(src)} -Ou -- -t MAF 2>/dev/null "
            f"| bcftools query -f '%MAF\n'", shell=True, executable="/bin/bash",
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        vals = [float(x) for x in proc.stdout.split() if x not in (".", "")]
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    n = len(vals)
    n_samples = len(samples) if samples else len(_vcf_sample_list(path))
    out = []
    ordered = sorted(marks)
    for i, m in enumerate(ordered):
        at = sum(1 for v in vals if v >= m)
        nxt = ordered[i + 1] if i + 1 < len(ordered) else None
        band = at - sum(1 for v in vals if v >= nxt) if nxt is not None else at
        out.append({"group": label, "n_samples": n_samples, "n_records": n,
                    "maf_min": m, "at_or_above": at,
                    "frac_at_or_above": round(at / n, 4) if n else 0.0,
                    "in_band_below_next": band})
    return out


def _vcf_sample_list(path: str) -> list[str]:
    return subprocess.run(["bcftools", "query", "-l", str(path)], stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, text=True).stdout.split()


def write_maf_spectrum(rows: list[dict], path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(MAF_SPECTRUM_COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in MAF_SPECTRUM_COLUMNS])


def maf_spectrum_note(rows: list[dict]) -> str:
    """One line for the whole callset, then one per group."""
    parts = []
    for label in dict.fromkeys(r["group"] for r in rows):
        rs = [r for r in rows if r["group"] == label]
        n = rs[0]["n_records"]
        marks = ", ".join(f"{r['maf_min']:.0%} {r['at_or_above']:,} ({r['frac_at_or_above']:.0%})"
                          for r in rs)
        parts.append(f"{label} (n={rs[0]['n_samples']}, {n:,} records): {marks}")
    return "; ".join(parts)


def write_variant_summary(rows: list[dict], path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(VARIANT_SUMMARY_COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in VARIANT_SUMMARY_COLUMNS])


def variant_summary_note(rows: list[dict]) -> str:
    """One line per class: ``snps 617 (biallelic 600 97.2%, triallelic 17 2.8%)``."""
    total = next((r["count"] for r in rows if r["class"] == "total"), 0)
    parts = []
    for r in rows:
        if r["class"] == "total" or r["n_alt"] != "all":
            continue
        subs = [f"{s['alleles']} {s['count']:,} {s['frac_class']:.1%}"
                for s in rows if s["class"] == r["class"] and s["n_alt"] != "all"]
        parts.append(f"{r['class']} {r['count']:,} {r['frac_total']:.1%} ({', '.join(subs)})")
    return f"{total:,} record(s)" + (": " + "; ".join(parts) if parts else "")


def sample_summary_table(path: str, *, ads_min: int = 10, frac_min: float = 0.80,
                         fws_min: float = 0.95, estimator: str = "regression",
                         min_depth: int = 0, n_bins: int = 10, min_alt_samples: int = 0,
                         snps_only: bool = True, multiallelic: str = "collapse",
                         trim: bool = True, exclude_call_regions=None) -> tuple[list[dict], int]:
    """Per-sample coverage and Fws, side by side, for a callset as it stands.

    The coverage columns are :func:`~.vcf_filters.sample_coverage_table`'s (the fraction of
    loci at ``ADS >= ads_min``, the margin against ``frac_min``) and the Fws columns are
    :func:`~.fws.fws_table`'s (the score, the sites it was scored over, whether it clears
    ``fws_min``). Nothing is dropped: ``would_drop_coverage`` and ``would_drop_fws`` say what
    the two filters *would* do at these thresholds, which is the point of running this on a
    final callset those filters may never have seen.

    ``FORMAT/ADS`` is what coverage is measured on; a callset without it (the filter that
    adds it switched off, or a file from elsewhere) gets the tag added on a temporary copy
    rather than being refused, since this changes nothing. Returns ``(rows, n_fws_sites)``.
    """
    from .fws import fws_table
    from .vcf_filters import _has_format_tag, sample_coverage_table

    require("bcftools")
    tmp = None
    src = path
    try:
        if not _has_format_tag(path, "ADS"):
            tmp = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
            sh(f"bcftools +fill-tags {q(path)} -Ob -o {q(tmp)} -- -t "
               f"{q('FORMAT/ADS=int(smpl_sum(FORMAT/AD))')}", tools=("bcftools",))
            src = tmp
        cov = {r["sample"]: r for r in sample_coverage_table(src, ads_min=ads_min,
                                                              frac_min=frac_min)}
        fws_rows, n_sites = fws_table(src, fws_min=fws_min, estimator=estimator,
                                      min_depth=min_depth, n_bins=n_bins,
                                      min_alt_samples=min_alt_samples, snps_only=snps_only,
                                      multiallelic=multiallelic, trim=trim,
                                      exclude_call_regions=exclude_call_regions)
    finally:
        if tmp and os.path.exists(tmp):
            os.unlink(tmp)
    fws = {r["sample"]: r for r in fws_rows}
    rows = []
    for name in cov:
        c, f = cov[name], fws.get(name, {})
        rows.append({
            "sample": name,
            "n_loci": c["n_loci"], "n_covered": c["n_covered"],
            "frac_covered": c["frac_covered"], "mean_ads": c["mean_ads"],
            "n_missing_ads": c["n_missing_ads"],
            "would_drop_coverage": c["dropped"],
            "fws": f.get("fws"), "fws_n_sites": f.get("n_sites", 0),
            "monoclonal": f.get("monoclonal", False),
            "would_drop_fws": f.get("dropped", True),
        })
    return sorted(rows, key=lambda r: (r["frac_covered"], r["sample"])), n_sites


SAMPLE_SUMMARY_COLUMNS = ["sample", "n_loci", "n_covered", "frac_covered", "mean_ads",
                          "n_missing_ads", "would_drop_coverage", "fws", "fws_n_sites",
                          "monoclonal", "would_drop_fws"]


def write_sample_summary(rows: list[dict], path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(SAMPLE_SUMMARY_COLUMNS)
        for r in rows:
            w.writerow(["" if r[c] is None else (f"{r[c]:.6f}" if c == "fws" else r[c])
                        for c in SAMPLE_SUMMARY_COLUMNS])


def sample_summary_note(rows: list[dict], n_sites: int, *, frac_min: float,
                        fws_min: float) -> str:
    n = len(rows)
    if not n:
        return "0 samples"
    thin = sum(r["would_drop_coverage"] for r in rows)
    scored = [r for r in rows if r["fws"] is not None]
    mono = sum(r["monoclonal"] for r in rows)
    med_cov = sorted(r["frac_covered"] for r in rows)[n // 2]
    med_fws = (sorted(r["fws"] for r in scored)[len(scored) // 2]) if scored else None
    return (f"{n} sample(s): median fraction covered {med_cov:.3f}, {thin} under {frac_min:g}; "
            f"Fws scored for {len(scored)} over {n_sites:,} site(s)"
            + (f", median {med_fws:.3f}" if med_fws is not None else "")
            + f", {mono} monoclonal at >= {fws_min:g}")
