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
"""

from __future__ import annotations

import json
from pathlib import Path

from . import vcf_filters as F
from .assets import resolve_bed
from .bcftools import count_variants
from .regenotype import filter_ad_regenotype


def _region(func):
    """Wrap a region filter so a ``builtin:`` bed value resolves to a shipped asset."""
    def run(inp, out, *, bed, **kw):
        return func(inp, out, bed=resolve_bed(bed), **kw)
    return run


# name -> callable(input_path, output_path, **params)
STEPS = {
    "hard_qc_filter": F.hard_qc_filter,
    "singleton_filter_add_ads": F.singleton_add_ads,
    "tandem_repeat_mask": _region(F.tandem_repeat_mask),
    "core_region_filter": _region(F.core_region_filter),
    "paralog_mask": _region(F.paralog_mask),
    "filter_ad_regenotype": lambda inp, out, **kw: filter_ad_regenotype(inp, out, **kw),
    "biallelic_snp_filter": F.biallelic_snp_filter,
    "sample_coverage_filter": F.sample_coverage_filter,
    "locus_missingness_filter": F.locus_missingness_filter,
    "maf_filter": F.maf_filter,
}

DEFAULT_CONFIG = {
    "steps": [
        {"name": "hard_qc_filter"},
        {"name": "singleton_filter_add_ads"},
        {"name": "tandem_repeat_mask", "params": {"bed": "builtin:pf3d7_tandem_repeats"}},
        {"name": "core_region_filter", "params": {"bed": "builtin:pf3d7_core_regions"}},
        {"name": "paralog_mask", "params": {"bed": "builtin:pf3d7_paralog_genes"}},
        {"name": "filter_ad_regenotype"},
        {"name": "biallelic_snp_filter"},
        {"name": "sample_coverage_filter"},
        {"name": "locus_missingness_filter"},
        {"name": "maf_filter", "params": {"maf_min": 0.02, "maf_max": 0.98}},
    ]
}


def load_config(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def run_pipeline(input_path: str, outdir: str, config: dict) -> list[dict]:
    """Run every step in order; return a per-step tally list."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    tally = [{"step": "input", "path": input_path, "variants": count_variants(input_path)}]
    prev = input_path
    for i, step in enumerate(config["steps"], start=1):
        name = step["name"]
        if name not in STEPS:
            raise SystemExit(f"ERROR: unknown pipeline step '{name}'. "
                             f"Known: {', '.join(STEPS)}")
        ext = step.get("ext", "bcf")
        params = step.get("params", {})
        out_path = str(out / f"{i:02d}_{name}.{ext}")
        print(f"[{i:02d}] {name} -> {out_path}")
        STEPS[name](prev, out_path, **params)
        n = count_variants(out_path)
        tally.append({"step": name, "path": out_path, "variants": n})
        print(f"     variants: {n:,}")
        prev = out_path

    return tally
