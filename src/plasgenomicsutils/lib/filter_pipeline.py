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
    "strip_stale_format": lambda inp, out, **kw: strip_stale_format(inp, out, **kw),
    "biallelic_snp_filter": F.biallelic_snp_filter,
    "sample_coverage_filter": F.sample_coverage_filter,
    "locus_missingness_filter": F.locus_missingness_filter,
    "maf_filter": F.maf_filter,
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
    "steps": [
        {"name": "hard_qc_filter"},
        # before singleton_filter_add_ads, which removes what this counts
        {"name": "singleton_counts", "report": True, "ext": "tsv",
         "params": {"min_depth": 5, "max_missing_frac": 0.2}},
        {"name": "singleton_filter_add_ads"},
        {"name": "tandem_repeat_mask", "params": {"bed": "builtin:pf3d7_tandem_repeats"}},
        {"name": "core_region_filter", "params": {"bed": "builtin:pf3d7_core_regions"}},
        {"name": "paralog_mask", "params": {"bed": "builtin:pf3d7_paralog_genes"}},
        {"name": "filter_ad_regenotype"},
        {"name": "biallelic_snp_filter"},
        {"name": "sample_coverage_filter"},
        {"name": "locus_missingness_filter"},
        {"name": "maf_filter", "params": {"maf_min": 0.02}},  # maf_max defaults to 1 - maf_min
    ]
}


def load_config(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def run_pipeline(input_path: str, outdir: str, config: dict,
                 *, emit_snp_bed: bool = True) -> list[dict]:
    """Run every step in order; return a per-step tally list.

    When ``emit_snp_bed`` (default), a BED of the final callset's SNPs is written next to
    the last step's output — the SNP panel the IBD tools read
    (``build_ibd_matrix --snp-format bed``).
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    tally = [{"step": "input", "path": input_path, "variants": count_variants(input_path)}]
    prev = input_path
    seen: list[str] = []
    for i, step in enumerate(config["steps"], start=1):
        name = step["name"]
        ext = step.get("ext", "bcf")
        params = step.get("params", {})
        out_path = str(out / f"{i:02d}_{name}.{ext}")

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
        STEPS[name](prev, out_path, **params)
        seen.append(name)
        index_vcf(out_path)   # keep intermediates indexed (quiets pysam, enables region queries)
        n = count_variants(out_path)
        tally.append({"step": name, "path": out_path, "variants": n})
        print(f"     variants: {n:,}")
        prev = out_path

    if emit_snp_bed and len(tally) > 1:
        bed_path = str(out / (Path(prev).stem + ".snps.bed"))
        F.snp_bed(prev, bed_path)
        n_snps = sum(1 for _ in open(bed_path))
        print(f"\nSNP panel BED ({n_snps:,} SNPs) -> {bed_path}")
        tally.append({"step": "snp_bed", "path": bed_path, "variants": n_snps})

    return tally
