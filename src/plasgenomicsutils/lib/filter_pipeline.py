"""Run an ordered, config-driven chain of filtering steps, tallying counts.

A pipeline config is JSON::

    {"steps": [
        {"name": "hard_qc_filter",         "params": {"qd": 20, "mq": 55}},
        {"name": "singleton_filter_add_ads"},
        {"name": "tandem_repeat_mask",     "params": {"bed": "tandems.bed"}, "ext": "vcf.gz"},
        {"name": "filter_ad_regenotype",   "params": {"min_reads": 2, "min_freq": 0.01}},
        {"name": "sample_coverage_filter"},
        {"name": "locus_missingness_filter"},
        {"name": "maf_filter",             "params": {"maf_min": 0.02, "maf_max": 0.98}}
    ]}

Each step writes ``<outdir>/NN_<name>.<ext>`` (ext defaults to ``bcf``) and its
input is the previous step's output.

A step marked ``"report": true`` is a **diagnostic**: it reads the current callset, writes
a table beside the filtered files, and passes the callset through untouched, so it can sit
anywhere in the chain without changing the result. Where it sits still matters -- see
:data:`REPORTS`.
"""

from __future__ import annotations

import json
from pathlib import Path

from . import vcf_filters as F
from .assets import resolve_bed
from .bcftools import count_variants, index_vcf
from .regenotype import filter_ad_regenotype
from .strip_format import strip_stale_format


def _region(func):
    """Wrap a region filter so a ``builtin:`` bed value resolves to a shipped asset.

    ``keep_bed`` (the whitelist) resolves the same way, though in practice it is a path to
    something dataset-specific rather than a shipped asset.
    """
    def run(inp, out, *, bed, keep_bed=None, **kw):
        return func(inp, out, bed=resolve_bed(bed),
                    keep_bed=resolve_bed(keep_bed) if keep_bed else None, **kw)
    return run


def _whitelisted(func, name):
    """Wrap a non-region step so its ``keep_bed`` resolves a ``builtin:`` value too."""
    def run(inp, out, *, keep_bed=None, **kw):
        return func(inp, out, keep_bed=resolve_bed(keep_bed) if keep_bed else None, **kw)
    return run


def _sidecar(out: str, name: str) -> str:
    """``09_sample_coverage_filter.bcf`` -> ``09_sample_coverage_filter_<name>``.

    A step's extra outputs sit next to its callset with the same numbered prefix, so the run
    directory reads in order and it is obvious which step produced what.
    """
    p = Path(out)
    stem = p.name
    for ext in (".vcf.gz", ".bcf.gz", ".vcf", ".bcf"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    return str(p.with_name(f"{stem}_{name}"))


def _sample_coverage(inp, out, **kw):
    """Drop low-coverage samples, leaving the coverage table that explains each drop.

    The table is written by default rather than on request: a sample vanishing from a cohort
    is exactly the kind of thing noticed weeks later, when re-running the step to find out why
    is expensive.
    """
    kw.setdefault("cov_table_path", _sidecar(out, "cov_info.tsv"))
    dropped = F.sample_coverage_filter(inp, out, **kw)
    print(f"     dropped {len(dropped)} low-coverage sample(s)"
          + (f": {', '.join(dropped)}" if dropped else "")
          + f"\n     coverage table -> {kw['cov_table_path']}")
    return dropped


# name -> callable(input_path, output_path, **params)
STEPS = {
    "no_alt_filter": _whitelisted(F.no_alt_filter, "no_alt_filter"),
    "hard_qc_filter": _whitelisted(F.hard_qc_filter, "hard_qc_filter"),
    "singleton_filter_add_ads": _whitelisted(F.singleton_add_ads, "singleton_filter_add_ads"),
    "tandem_repeat_mask": _region(F.tandem_repeat_mask),
    "core_region_filter": _region(F.core_region_filter),
    "paralog_mask": _region(F.paralog_mask),
    "filter_ad_regenotype": lambda inp, out, **kw: filter_ad_regenotype(inp, out, **kw),
    "strip_stale_format": lambda inp, out, **kw: strip_stale_format(inp, out, **kw),
    "biallelic_snp_filter": F.biallelic_snp_filter,
    "sample_coverage_filter": lambda inp, out, **kw: _sample_coverage(inp, out, **kw),
    "locus_missingness_filter": _whitelisted(F.locus_missingness_filter,
                                             "locus_missingness_filter"),
    "maf_filter": _whitelisted(F.maf_filter, "maf_filter"),
}

def _singleton_report(inp, out, **kw):
    """Per-sample singleton counts + near-identical pairs, written to `out`."""
    from .singletons import count_singletons, flag_outliers

    mad = kw.pop("mad_cutoff", None)
    dup = kw.pop("duplicate_frac", None)
    df, n_variants = count_singletons(inp, **kw)
    df = flag_outliers(df, **{k: v for k, v in
                              (("mad_cutoff", mad), ("duplicate_frac", dup))
                              if v is not None})
    df = df.sort_values("singleton_rate", ascending=False)
    df.to_csv(out, sep="\t", index=False)
    flagged = df[df["flag"] != ""]
    print(f"     {n_variants:,} variants scanned, "
          f"{int(df['n_singleton'].sum()):,} singletons, {len(flagged)} sample(s) flagged")
    for r in flagged.itertuples(index=False):
        print(f"       {r.sample}\t{r.singleton_rate:.2f}/1000\t{r.flag}")
    return len(df)


#: Steps that take a ``keep_bed`` whitelist. The rest either judge samples rather than
#: variants (``sample_coverage_filter``), transform genotypes without dropping records
#: (``filter_ad_regenotype``, ``strip_stale_format``), or define what the callset *is*
#: rather than filtering it on quality (``biallelic_snp_filter`` -- letting a whitelisted
#: multiallelic record through would break every downstream reader's assumption).
WHITELISTABLE = {
    "no_alt_filter", "hard_qc_filter", "singleton_filter_add_ads", "tandem_repeat_mask",
    "core_region_filter", "paralog_mask", "locus_missingness_filter", "maf_filter",
}

#: name -> callable(input_path, output_path, **params). A report reads the callset and
#: writes a table; it never changes the data.
#:
#: ``singleton_counts`` has to run **before** ``singleton_filter_add_ads``, which drops
#: exactly the variants it counts -- run it after and every sample scores zero. The
#: default config places it right after ``hard_qc_filter``, so obvious junk is gone but
#: the private variants are still there.
REPORTS = {
    "singleton_counts": _singleton_report,
}

DEFAULT_CONFIG = {
    # A whitelist BED applied to every step in WHITELISTABLE: the variants it covers survive
    # those filters. Written out set to null so the key is discoverable -- JSON has no
    # comments, and a resistance locus that a MAF floor or a coverage rule would otherwise
    # remove is exactly the thing worth keeping. A step's own params.keep_bed overrides it.
    "keep_bed": None,
    "steps": [
        # Non-variant records first, and in their own step: the bias statistics are
        # computed whether or not an ALT was called, so hard_qc_filter does remove them --
        # counting them separately keeps "nothing to call here" apart from "failed QC".
        {"name": "no_alt_filter", "params": {"keep": False}},
        # `caller` is written out at its default for the same reason as keep_bed: a
        # bcftools callset carries none of GATK's metrics, and "bcftools" here is what
        # makes this step read the ones it does carry (FS/RPBZ/SCBZ/MQBZ/MQSBZ).
        {"name": "hard_qc_filter", "params": {"caller": "gatk"}},
        # before singleton_filter_add_ads, which removes what this counts
        {"name": "singleton_counts", "report": True, "ext": "tsv",
         "params": {"min_depth": 5, "max_missing_frac": 0.2}},
        {"name": "singleton_filter_add_ads"},
        {"name": "tandem_repeat_mask", "params": {"bed": "builtin:pf3d7_tandem_repeats"}},
        {"name": "core_region_filter", "params": {"bed": "builtin:pf3d7_core_regions"}},
        # Off by default. core_region_filter has already dropped the subtelomeric multigene
        # families that mismap worst; most of what this removes next sits in the core and is
        # single-copy, so some of it genuinely misbehaves but a lot of it is fine. Losing all
        # of it costs real signal, and paralogy is not by itself evidence a locus is unusable.
        # Set "enabled": true when mismapping is the thing being controlled for.
        {"name": "paralog_mask", "enabled": False,
         "params": {"bed": "builtin:pf3d7_paralog_genes"}},
        {"name": "filter_ad_regenotype"},
        {"name": "biallelic_snp_filter"},
        {"name": "sample_coverage_filter"},
        {"name": "locus_missingness_filter"},
        {"name": "maf_filter", "params": {"maf_min": 0.02}},  # maf_max defaults to 1 - maf_min
    ]
}


def load_config(path: str) -> dict:
    """Read a pipeline config from JSON."""
    with open(path) as fh:
        return json.load(fh)


def run_pipeline(input_path: str, outdir: str, config: dict,
                 *, emit_snp_bed: bool = True) -> list[dict]:
    """Run every step in order; return a per-step tally list.

    When ``emit_snp_bed`` (default), a BED of the final callset's SNPs is written next to
    the last step's output — the SNP panel the IBD tools read
    (``build_ibd_matrix --snp-format bed``).

    A top-level ``"keep_bed"`` in the config is a whitelist for the whole chain: every step
    in :data:`WHITELISTABLE` keeps the variants it covers, so a resistance locus survives the
    run without being written into each step. A step's own ``params.keep_bed`` overrides it,
    and ``"keep_bed": null`` in a step's params opts that step out. A whitelist that rescues
    nothing is warned about once for the run rather than at each step, since a step with
    nothing to rescue is the normal case and most steps are that.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    tally = [{"step": "input", "path": input_path, "variants": count_variants(input_path)}]
    prev = input_path
    seen: list[str] = []
    # One "the whitelist rescued nothing" warning for the run, not one per step: most
    # steps have nothing for a whitelist to do, and saying so each time reads as an error.
    with F.deferred_whitelist_warnings():
        for i, step in enumerate(config["steps"], start=1):
            name = step["name"]
            ext = step.get("ext", "bcf")
            params = dict(step.get("params", {}))
            if (config.get("keep_bed") and name in WHITELISTABLE
                    and "keep_bed" not in params):
                params["keep_bed"] = config["keep_bed"]
            out_path = str(out / f"{i:02d}_{name}.{ext}")

            # `"enabled": false` keeps a step in the config, and out of the run. That is how an
            # optional step stays discoverable -- JSON has no comments, so a step you would
            # otherwise have to know about is written down with the switch off.
            if step.get("enabled", True) is False:
                print(f"[{i:02d}] {name} -- skipped (\"enabled\": false)")
                tally.append({"step": name, "skipped": True})
                continue

            if step.get("report"):
                if name not in REPORTS:
                    raise SystemExit(f"ERROR: unknown pipeline report '{name}'. "
                                     f"Known: {', '.join(REPORTS)}")
                if name == "singleton_counts" and "singleton_filter_add_ads" in seen:
                    print(f"[{i:02d}] WARNING: singleton_counts runs after "
                          f"singleton_filter_add_ads, which drops the variants it counts -- "
                          f"every sample will score zero. Move it earlier.")
                print(f"[{i:02d}] {name} (report) -> {out_path}")
                n = REPORTS[name](prev, out_path, **params)
                tally.append({"step": name, "path": out_path, "report": True, "rows": n})
                seen.append(name)
                continue                       # the callset is unchanged

            if name not in STEPS:
                raise SystemExit(f"ERROR: unknown pipeline step '{name}'. "
                                 f"Known: {', '.join(STEPS)}")
            print(f"[{i:02d}] {name} -> {out_path}")
            rescued = STEPS[name](prev, out_path, **params)
            seen.append(name)
            index_vcf(out_path)   # keep intermediates indexed (quiets pysam, enables region queries)
            n = count_variants(out_path)
            row = {"step": name, "path": out_path, "variants": n}
            if isinstance(rescued, int) and rescued > 0:
                row["rescued"] = rescued
            tally.append(row)
            print(f"     variants: {n:,}")
            prev = out_path

    if emit_snp_bed and len(tally) > 1:
        bed_path = str(out / (Path(prev).stem + ".snps.bed"))
        F.snp_bed(prev, bed_path)
        n_snps = sum(1 for _ in open(bed_path))
        print(f"\nSNP panel BED ({n_snps:,} SNPs) -> {bed_path}")
        tally.append({"step": "snp_bed", "path": bed_path, "variants": n_snps})

    return tally
