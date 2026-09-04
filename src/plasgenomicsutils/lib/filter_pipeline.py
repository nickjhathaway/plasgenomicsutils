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

import difflib
import importlib
import os
import inspect
import json
from datetime import datetime
from pathlib import Path

from . import vcf_filters as F
from .assets import resolve_bed
from .. import __version__
from .bcftools import VARIANT_TYPES, index_vcf, variant_type_counts
from .regenotype import filter_ad_regenotype
from .reporting import detail, listing, say
from .strip_format import strip_stale_format


def _region(func):
    """Wrap a region filter so a ``builtin:`` bed value resolves to a shipped asset.

    ``keep_bed`` (the whitelist) resolves the same way, though in practice it is a path to
    something dataset-specific rather than a shipped asset.
    """
    def run(inp, out, *, bed, keep_bed=None, **kw):
        return func(inp, out, bed=resolve_bed(bed),
                    keep_bed=resolve_bed(keep_bed) if keep_bed else None, **kw)
    run.target = func            # what `params` is validated against
    return run


def _whitelisted(func, name):
    """Wrap a non-region step so its ``keep_bed`` resolves a ``builtin:`` value too."""
    def run(inp, out, *, keep_bed=None, **kw):
        return func(inp, out, keep_bed=resolve_bed(keep_bed) if keep_bed else None, **kw)
    run.target = func            # what `params` is validated against
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
    say(f"     dropped {len(dropped)} low-coverage sample(s)" + listing(dropped, limit=5)
        + f"\n     coverage table -> {kw['cov_table_path']}")
    return dropped


def _fws(inp, out, **kw):
    """Keep the monoclonal samples, leaving the Fws table that explains each drop."""
    from .fws import fws_filter

    kw.setdefault("fws_table_path", _sidecar(out, "fws.tsv"))
    dropped = fws_filter(inp, out, **kw)
    say(f"     dropped {len(dropped)} polyclonal/unscored sample(s)"
        + listing(dropped, limit=5) + f"\n     Fws table -> {kw['fws_table_path']}")
    return dropped


_sample_coverage.target = F.sample_coverage_filter
_fws.target_ref = (".fws", "fws_filter")


# name -> callable(input_path, output_path, **params)
STEPS = {
    "no_alt_filter": _whitelisted(F.no_alt_filter, "no_alt_filter"),
    "hard_qc_filter": _whitelisted(F.hard_qc_filter, "hard_qc_filter"),
    "singleton_filter_add_ads": _whitelisted(F.singleton_add_ads, "singleton_filter_add_ads"),
    "tandem_repeat_mask": _region(F.tandem_repeat_mask),
    "core_region_filter": _region(F.core_region_filter),
    "paralog_mask": _region(F.paralog_mask),
    "filter_ad_regenotype": filter_ad_regenotype,
    "strip_stale_format": strip_stale_format,
    "biallelic_snp_filter": F.biallelic_snp_filter,
    "sample_coverage_filter": _sample_coverage,
    "fws_filter": _fws,
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
    say(f"     {n_variants:,} variants scanned, "
        f"{int(df['n_singleton'].sum()):,} singletons, {len(flagged)} sample(s) flagged")
    for r in flagged.itertuples(index=False):
        detail(f"       {r.sample}\t{r.singleton_rate:.2f}/1000\t{r.flag}")
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
_singleton_report.target_ref = (".singletons", "count_singletons")
_singleton_report.extra_params = ("mad_cutoff", "duplicate_frac")

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
        # trim first, so every count below describes the alleles this cohort actually
        # carries rather than the ones the callset it was subset from did
        {"name": "no_alt_filter", "params": {"keep": False, "trim": True}},
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
        # both tests written out at their defaults so the split is discoverable: `biallelic`
        # off keeps multiallelic SNPs for downstream tools that can read them.
        {"name": "biallelic_snp_filter",
         "params": {"snps_only": True, "biallelic": True, "mnp_handling": "split"}},
        {"name": "sample_coverage_filter"},
        {"name": "locus_missingness_filter"},
        {"name": "maf_filter", "params": {"maf_min": 0.02}},  # maf_max defaults to 1 - maf_min
        # Keeping only monoclonal infections is an analysis choice, not a QC rule -- it
        # changes which infections the callset describes -- so it is written out switched
        # off rather than left undiscoverable. It drops samples and no variants; a site the
        # survivors no longer support is still there, which is why re-running maf_filter and
        # locus_missingness_filter after it is worth doing when the frequencies matter.
        {"name": "fws_filter", "enabled": False, "params": {"fws_min": 0.95}},
    ]
}


def load_config(path: str) -> dict:
    """Read a pipeline config from JSON."""
    with open(path) as fh:
        return json.load(fh)


#: Keys that belong on a step rather than inside its ``params``. Putting one in ``params``
#: is a slip the step cannot catch on its own -- it just arrives as an unexpected argument.
STEP_KEYS = ("name", "params", "ext", "enabled", "report")

#: Keys a config may carry at the top level.
CONFIG_KEYS = ("steps", "keep_bed", "remove_intermediates", "_meta")


def _accepted_params(step) -> set[str] | None:
    """The params a step takes, or ``None`` when it takes anything (``**kwargs``).

    The first two parameters are the input and output paths, whatever they are called, so
    they are dropped by position rather than by name.
    """
    fn = getattr(step, "target", None)
    if fn is None:
        # a step whose real work is imported inside its body names it instead, so the import
        # stays where it was put rather than being hoisted for the sake of validation
        ref = getattr(step, "target_ref", None)
        fn = (getattr(importlib.import_module(ref[0], __package__), ref[1]) if ref else step)
    try:
        params = list(inspect.signature(fn).parameters.values())
    except (TypeError, ValueError):          # pragma: no cover - builtins, C callables
        return None
    if any(p.kind is p.VAR_KEYWORD for p in params):
        return None
    kinds = (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY)
    return ({p.name for p in params[2:] if p.kind in kinds}
            | set(getattr(step, "extra_params", ())))


def validate_config(config: dict) -> None:
    """Check a config before anything runs.

    A pipeline writes files as it goes, so a mistake in step nine is worth catching before
    step one rather than after eight steps of output. Everything here is knowable from the
    config alone.
    """
    unknown = sorted(set(config) - set(CONFIG_KEYS))
    if unknown:
        near = difflib.get_close_matches(unknown[0], CONFIG_KEYS, n=1)
        raise SystemExit(
            f"ERROR: unknown top-level config key '{unknown[0]}'"
            + (f". Did you mean '{near[0]}'?" if near else "")
            + f"\n  The config takes: {', '.join(CONFIG_KEYS)}")
    steps = config.get("steps")
    if not isinstance(steps, list) or not steps:
        raise SystemExit("ERROR: the config needs a non-empty \"steps\" list")
    for i, step in enumerate(steps, start=1):
        where = f"step {i}"
        if not isinstance(step, dict) or "name" not in step:
            raise SystemExit(f"ERROR: {where} has no \"name\"")
        name = step["name"]
        where = f"{where} ({name})"
        table = REPORTS if step.get("report") else STEPS
        if name not in table:
            kind = "report" if step.get("report") else "step"
            near = difflib.get_close_matches(name, table, n=1)
            raise SystemExit(
                f"ERROR: {where}: unknown pipeline {kind} '{name}'"
                + (f". Did you mean '{near[0]}'?" if near else "")
                + f"\n  Known: {', '.join(sorted(table))}")

        params = step.get("params") or {}
        if not isinstance(params, dict):
            raise SystemExit(f"ERROR: {where}: \"params\" must be an object")
        accepted = _accepted_params(table[name])
        for key in params:
            if key in STEP_KEYS:
                raise SystemExit(
                    f"ERROR: {where}: \"{key}\" goes on the step, not inside its params -- "
                    f"it is a pipeline control, not an argument to {name}.\n"
                    f'  {{"name": "{name}", "{key}": ..., "params": {{...}}}}')
            if accepted is not None and key not in accepted:
                near = difflib.get_close_matches(key, accepted, n=1)
                raise SystemExit(
                    f"ERROR: {where}: unknown param '{key}'"
                    + (f". Did you mean '{near[0]}'?" if near else "")
                    + (f"\n  {name} takes: {', '.join(sorted(accepted))}" if accepted
                       else f"\n  {name} takes no params"))


def _jsonable(v):
    """A default rendered for JSON, so a config record is a config you can run again."""
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, (set, frozenset)):
        return sorted(_jsonable(x) for x in v)
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in v.items()}
    return str(v)


def _step_defaults(step) -> dict:
    """Every parameter a step has a default for, read from its own signature."""
    fn = getattr(step, "target", None)
    if fn is None:
        ref = getattr(step, "target_ref", None)
        fn = (getattr(importlib.import_module(ref[0], __package__), ref[1]) if ref else step)
    try:
        params = list(inspect.signature(fn).parameters.values())
    except (TypeError, ValueError):          # pragma: no cover - builtins, C callables
        return {}
    kinds = (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY)
    return {p.name: _jsonable(p.default) for p in params[2:]
            if p.kind in kinds and p.default is not inspect.Parameter.empty}


def effective_config(config: dict, **meta) -> dict:
    """The config as it will actually run, with every unset parameter filled in.

    A config only records what somebody chose to write down, so the thresholds that did the
    work are mostly absent from it -- they live in the code. Resolving them into the run's own
    directory means the answer to "what cutoffs were these?" is a file beside the results
    rather than a version of the package you have to go and find.

    Steps switched off keep their entry and their defaults: what did **not** run is part of
    the record too. The result is a valid config -- run it again and you get this run.
    """
    out = {k: config[k] for k in ("keep_bed", "remove_intermediates") if k in config}
    out.setdefault("keep_bed", None)
    out.setdefault("remove_intermediates", False)
    if meta:
        out["_meta"] = {k: _jsonable(v) for k, v in meta.items()}
    steps = []
    for step in config["steps"]:
        name = step["name"]
        table = REPORTS if step.get("report") else STEPS
        params = dict(_step_defaults(table[name]))
        params.pop("keep_bed", None)
        if name in WHITELISTABLE:
            params["keep_bed"] = config.get("keep_bed")
        params.update(step.get("params") or {})
        entry = {"name": name}
        for key in ("report", "ext", "enabled"):
            if key in step:
                entry[key] = step[key]
        entry["params"] = params
        steps.append(entry)
    out["steps"] = steps
    return out


def type_counts_note(counts: dict) -> str:
    """Every class present, SNPs first. Used by the per-step lines and the end-of-run
    table alike, so the two cannot drift into showing different things."""
    if not counts:
        return ""
    named = [f"snps {counts.get('snps', 0):,}"]
    named += [f"{n} {counts[n]:,}" for n in VARIANT_TYPES if n != "snps" and counts.get(n)]
    return "   (" + ", ".join(named) + ")"


def _types_note(counts: dict) -> str:
    return type_counts_note(counts)


def _remove_intermediate(row: dict) -> None:
    """Delete a step's callset and index once the next step has read it."""
    path = row.get("path")
    if not path:
        return
    for suffix in ("", ".csi", ".tbi"):
        try:
            os.remove(path + suffix)
        except OSError:
            pass
    row["removed"] = True
    say(f"     removed {Path(path).name}")


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
    validate_config(config)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    used = out / "config_used.json"
    used.write_text(json.dumps(
        effective_config(config, version=__version__, input=os.path.abspath(input_path),
                         started=datetime.now().astimezone().isoformat(timespec="seconds")),
        indent=2) + "\n")
    say(f"  config as run -> {used}")

    counts = variant_type_counts(input_path)
    tally = [{"step": "input", "path": input_path, "variants": counts["total"],
              "types": counts}]
    prev = input_path
    seen: list[str] = []
    # `remove_intermediates` deletes each step's callset as soon as the next one has read it,
    # so a long chain over a large cohort costs one intermediate on disk rather than all of
    # them. The input is never touched, the final output stays, and the side tables -- the
    # record of what happened -- are small and are kept whatever this says.
    prune = bool(config.get("remove_intermediates", False))
    prev_row: dict | None = None
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
                say(f"[{i:02d}] {name} -- skipped (\"enabled\": false)")
                tally.append({"step": name, "skipped": True})
                continue

            if step.get("report"):
                if name not in REPORTS:
                    raise SystemExit(f"ERROR: unknown pipeline report '{name}'. "
                                     f"Known: {', '.join(REPORTS)}")
                if name == "singleton_counts" and "singleton_filter_add_ads" in seen:
                    say(f"[{i:02d}] WARNING: singleton_counts runs after "
                          f"singleton_filter_add_ads, which drops the variants it counts -- "
                          f"every sample will score zero. Move it earlier.")
                say(f"[{i:02d}] {name} (report) -> {out_path}")
                n = REPORTS[name](prev, out_path, **params)
                tally.append({"step": name, "path": out_path, "report": True, "rows": n})
                seen.append(name)
                continue                       # the callset is unchanged

            if name not in STEPS:
                raise SystemExit(f"ERROR: unknown pipeline step '{name}'. "
                                 f"Known: {', '.join(STEPS)}")
            say(f"[{i:02d}] {name} -> {out_path}")
            rescued = STEPS[name](prev, out_path, **params)
            seen.append(name)
            index_vcf(out_path)   # keep intermediates indexed (quiets pysam, enables region queries)
            counts = variant_type_counts(out_path)
            n = counts["total"]
            row = {"step": name, "path": out_path, "variants": n, "types": counts}
            if isinstance(rescued, int) and rescued > 0:
                row["rescued"] = rescued
            tally.append(row)
            say(f"     variants: {n:,}{_types_note(counts)}")
            if prune and prev_row is not None:
                _remove_intermediate(prev_row)
            prev_row = row
            prev = out_path

    if emit_snp_bed and len(tally) > 1:
        bed_path = str(out / (Path(prev).stem + ".snps.bed"))
        F.snp_bed(prev, bed_path)
        n_snps = sum(1 for _ in open(bed_path))
        say(f"\nSNP panel BED ({n_snps:,} SNPs) -> {bed_path}")
        tally.append({"step": "snp_bed", "path": bed_path, "variants": n_snps})

    return tally
