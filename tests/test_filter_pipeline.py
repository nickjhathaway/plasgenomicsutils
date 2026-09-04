"""The filter pipeline's report steps: diagnostics that read the callset but never
change it, and the ordering trap that makes one of them meaningless."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"
BCF = DATA / "ghana_cambodia.pf7.tiny.bcf"


@pytest.fixture
def with_ads(tmp_path):
    """The fixture callset with FORMAT/ADS added -- sample coverage is measured on ADS, so
    the filter needs the tag that `singleton_filter_add_ads` writes."""
    from plasgenomicsutils.lib import vcf_filters as F

    out = tmp_path / "with_ads.bcf"
    F.singleton_add_ads(str(BCF), str(out))
    return str(out)


# --------------------------------------------------------------------------- #
#  Report steps                                                                #
# --------------------------------------------------------------------------- #


def test_a_report_step_writes_a_table_and_leaves_the_callset_alone(tmp_path):
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    cfg = {"steps": [
        {"name": "singleton_counts", "report": True, "ext": "tsv",
         "params": {"min_depth": 0, "max_missing_frac": 0.5}},
        {"name": "singleton_filter_add_ads"},
    ]}
    tally = run_pipeline(str(BCF), str(tmp_path), cfg, emit_snp_bed=False)
    report = next(t for t in tally if t.get("report"))
    assert Path(report["path"]).exists()
    assert report["rows"] == 60                       # one row per sample
    # the filter step read the original input, not the report's path
    filt = next(t for t in tally if t["step"] == "singleton_filter_add_ads")
    assert filt["variants"] < tally[0]["variants"]


def test_counting_singletons_after_the_singleton_filter_is_warned_about(tmp_path, capsys):
    """The filter drops exactly what the report counts, so the order is a real trap."""
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    cfg = {"steps": [
        {"name": "singleton_filter_add_ads"},
        {"name": "singleton_counts", "report": True, "ext": "tsv"},
    ]}
    run_pipeline(str(BCF), str(tmp_path), cfg, emit_snp_bed=False)
    out = capsys.readouterr().out
    assert "WARNING" in out and "Move it earlier" in out
    assert "0 singletons" in out


def test_the_default_config_counts_singletons_before_filtering_them():
    from plasgenomicsutils.lib.filter_pipeline import DEFAULT_CONFIG

    names = [s["name"] for s in DEFAULT_CONFIG["steps"]]
    assert names.index("singleton_counts") < names.index("singleton_filter_add_ads")


def test_an_unknown_report_name_is_rejected(tmp_path):
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    cfg = {"steps": [{"name": "nope", "report": True}]}
    with pytest.raises(SystemExit, match="unknown pipeline report"):
        run_pipeline(str(BCF), str(tmp_path), cfg, emit_snp_bed=False)


# --------------------------------------------------------------------------- #
#  Optional steps                                                              #
# --------------------------------------------------------------------------- #


def test_a_step_marked_disabled_is_skipped_and_does_not_break_the_chain(tmp_path):
    """`"enabled": false` keeps a step visible in the config and out of the run."""
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    cfg = {"steps": [
        {"name": "paralog_mask", "enabled": False,
         "params": {"bed": "builtin:pf3d7_paralog_genes"}},
        {"name": "singleton_filter_add_ads"},
    ]}
    tally = run_pipeline(str(BCF), str(tmp_path), cfg, emit_snp_bed=False)
    skipped = next(t for t in tally if t["step"] == "paralog_mask")
    assert skipped["skipped"] is True
    assert "path" not in skipped                    # nothing was written
    # the next step read the original input, not a missing paralog_mask output
    filt = next(t for t in tally if t["step"] == "singleton_filter_add_ads")
    assert filt["variants"] < tally[0]["variants"]


def test_paralog_masking_is_off_by_default_but_present_to_turn_on():
    from plasgenomicsutils.lib.filter_pipeline import DEFAULT_CONFIG

    step = next(s for s in DEFAULT_CONFIG["steps"] if s["name"] == "paralog_mask")
    assert step["enabled"] is False              # discoverable, not silently absent
    assert step["params"]["bed"] == "builtin:pf3d7_paralog_genes"
    # exactly two steps ship switched off, and both for a stated reason: paralog masking is
    # a choice about which regions to trust, and fws_filter changes which infections the
    # callset describes. Anything else arriving here off is a mistake, not a default.
    off = {s["name"] for s in DEFAULT_CONFIG["steps"] if s.get("enabled", True) is False}
    assert off == {"paralog_mask", "fws_filter"}


def test_the_default_config_still_masks_tandem_repeats_and_keeps_the_core():
    """Turning paralog masking off must not touch the masks that stayed on."""
    from plasgenomicsutils.lib.filter_pipeline import DEFAULT_CONFIG

    on = [s["name"] for s in DEFAULT_CONFIG["steps"] if s.get("enabled", True)]
    assert "tandem_repeat_mask" in on and "core_region_filter" in on


# --------------------------------------------------------------------------- #
#  The per-step summary                                                        #
# --------------------------------------------------------------------------- #


def test_the_summary_handles_filter_report_and_skipped_rows():
    """A report row records `rows` and a skipped row records neither, so the summary
    cannot just read `variants` off every row."""
    from plasgenomicsutils.scripts.vcf.filter_pipeline import _tally_fields

    assert _tally_fields({"step": "a", "path": "p", "variants": 7}) == ("variants", 7, "p")
    assert _tally_fields({"step": "b", "path": "p", "report": True, "rows": 60}) == \
        ("report", 60, "p")
    assert _tally_fields({"step": "c", "skipped": True}) == ("skipped", "", "")


def test_the_default_pipeline_writes_a_summary_for_every_row(tmp_path):
    """End to end through the CLI wrapper: a report step and a disabled step both have to
    survive the summary, which is what reading `variants` off every row used to break."""
    import json
    import sys

    from plasgenomicsutils.scripts.vcf import filter_pipeline as FP

    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"steps": [
        {"name": "singleton_counts", "report": True, "ext": "tsv",
         "params": {"min_depth": 0, "max_missing_frac": 0.5}},
        {"name": "paralog_mask", "enabled": False,
         "params": {"bed": "builtin:pf3d7_paralog_genes"}},
        {"name": "singleton_filter_add_ads"},
    ]}))
    outdir = tmp_path / "out"
    argv = sys.argv
    sys.argv = ["filter_pipeline", "--input", str(BCF), "--config", str(cfg),
                "--outdir", str(outdir), "--no-snp-bed"]
    try:
        FP.filter_pipeline()
    finally:
        sys.argv = argv

    rows = (outdir / "variant_counts.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    assert header[:3] == ["step", "kind", "count"]
    assert header[-3:] == ["rescued", "removed", "path"]
    assert header[3:-3] == ["snps", "indels", "mnps", "mixed", "other", "no_alt"]
    cols = {n: i for i, n in enumerate(header)}
    # nothing is whitelisted in this config, and nothing was pruned
    assert all(r.split("\t")[cols["rescued"]] == "" for r in rows[1:])
    assert all(r.split("\t")[cols["removed"]] == "" for r in rows[1:])
    kinds = {r.split("\t")[0]: r.split("\t")[1] for r in rows[1:]}
    assert kinds["singleton_counts"] == "report"
    assert kinds["paralog_mask"] == "skipped"
    assert kinds["singleton_filter_add_ads"] == "variants"


# --------------------------------------------------------------------------- #
#  Sample-coverage diagnostics                                                 #
# --------------------------------------------------------------------------- #


def test_the_coverage_table_accounts_for_every_keep_and_drop(tmp_path, with_ads):
    """One row per sample, and `dropped` agrees with the frac_covered/frac_min comparison
    the filter makes -- the table has to explain the decision, not merely resemble it."""
    from plasgenomicsutils.lib import vcf_filters as F

    rows = F.sample_coverage_table(with_ads, ads_min=5, frac_min=0.5)
    assert len(rows) == 60
    assert [r["frac_covered"] for r in rows] == sorted(r["frac_covered"] for r in rows)
    for r in rows:
        assert r["dropped"] == (r["frac_covered"] < r["frac_min"])
        assert r["dropped"] == (r["margin"] < 0)
        assert r["n_covered"] <= r["n_loci"]
        assert r["ads_min"] == 5 and r["frac_min"] == 0.5

    out = tmp_path / "cov.bcf"
    dropped = F.sample_coverage_filter(with_ads, str(out), ads_min=5, frac_min=0.5)
    assert dropped == sorted(r["sample"] for r in rows if r["dropped"])


def test_a_threshold_nobody_can_meet_drops_everyone_and_says_so(with_ads):
    from plasgenomicsutils.lib import vcf_filters as F

    rows = F.sample_coverage_table(with_ads, ads_min=10_000, frac_min=0.5)
    assert all(r["n_covered"] == 0 and r["dropped"] for r in rows)
    # ... and nothing is dropped when the bar is on the floor
    rows = F.sample_coverage_table(with_ads, ads_min=1, frac_min=0.0)
    assert not any(r["dropped"] for r in rows)


def test_a_callset_without_ADS_is_refused_rather_than_emptied(tmp_path):
    """Coverage is measured on FORMAT/ADS. Without it every sample scores zero and the whole
    cohort would be dropped, which is much worse than stopping."""
    from plasgenomicsutils.lib import vcf_filters as F

    with pytest.raises(SystemExit, match="no FORMAT/ADS"):
        F.sample_coverage_filter(str(BCF), str(tmp_path / "o.bcf"))


def test_the_pipeline_writes_the_coverage_table_beside_the_step(tmp_path):
    """Named for the step that produced it, so the run directory reads in order."""
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    cfg = {"steps": [
        {"name": "singleton_filter_add_ads"},
        {"name": "sample_coverage_filter", "params": {"ads_min": 5, "frac_min": 0.5}},
    ]}
    run_pipeline(str(BCF), str(tmp_path), cfg, emit_snp_bed=False)
    tsv = tmp_path / "02_sample_coverage_filter_cov_info.tsv"
    assert tsv.exists()

    lines = tsv.read_text().splitlines()
    assert lines[0].split("\t") == [
        "sample", "n_loci", "n_covered", "frac_covered", "n_missing_ads", "mean_ads",
        "ads_min", "frac_min", "margin", "dropped"]
    assert len(lines) == 61                                   # header + 60 samples

    # the samples marked dropped are exactly the ones missing from the output callset
    from plasgenomicsutils.lib import vcf_filters as F
    marked = {r.split("\t")[0] for r in lines[1:] if r.split("\t")[-1] == "True"}
    remaining = set(F._samples(str(tmp_path / "02_sample_coverage_filter.bcf")))
    assert marked and marked.isdisjoint(remaining)
    assert marked | remaining == set(F._samples(str(BCF)))


def test_a_sidecar_name_survives_a_two_part_extension():
    from plasgenomicsutils.lib.filter_pipeline import _sidecar

    assert _sidecar("out/09_sample_coverage_filter.bcf", "cov_info.tsv") == \
        "out/09_sample_coverage_filter_cov_info.tsv"
    assert _sidecar("out/09_step.vcf.gz", "cov_info.tsv") == "out/09_step_cov_info.tsv"


def test_borderline_samples_are_printed_so_a_surprise_drop_is_visible(tmp_path, with_ads,
                                                                     capsys):
    """A sample that missed by a hair reads very differently from a hopeless one, so both
    sides of the threshold are logged, not just the drops."""
    from plasgenomicsutils.lib import vcf_filters as F

    F.sample_coverage_filter(with_ads, str(tmp_path / "o.bcf"), ads_min=5, frac_min=0.5)
    out = capsys.readouterr().out
    assert "DROPPED" in out and "kept" in out
    assert "mean ADS" in out                      # why, not just whether


# --------------------------------------------------------------------------- #
#  Region whitelists                                                           #
# --------------------------------------------------------------------------- #


@pytest.fixture
def region_case(tmp_path):
    """The first four fixture records, a mask covering all of them, and a one-base whitelist
    on the second and third. Those two have multi-base REFs, which is the case that matters."""
    import subprocess

    rows = subprocess.run(
        f"bcftools query -f '%CHROM\\t%POS\\t%REF\\n' {BCF} | head -4",
        shell=True, capture_output=True, text=True).stdout.strip().splitlines()
    recs = [(r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[2]) for r in rows]
    chrom, picks = recs[0][0], [p for _, p, _ in recs]
    assert any(len(ref) > 1 for _, _, ref in recs)      # the point of the fixture

    mask = tmp_path / "mask.bed"
    mask.write_text(f"{chrom}\t{min(picks) - 1}\t{max(picks) + 5}\n")
    keep = tmp_path / "keep.bed"
    keep.write_text("".join(f"{chrom}\t{p - 1}\t{p}\n" for p in picks[1:3]))
    return {"chrom": chrom, "picks": picks, "mask": str(mask), "keep": str(keep),
            "rescued": picks[1:3]}


def _positions(path, among):
    import subprocess

    out = subprocess.run(f"bcftools query -f '%POS\\n' {path}", shell=True,
                         capture_output=True, text=True).stdout.split()
    return [int(x) for x in out if int(x) in among]


def test_a_one_base_whitelist_rescues_a_multi_base_record_from_a_drop_mask(tmp_path, region_case):
    """bedtools sizes a record as [POS-1, POS-1+len(REF)), so carving only the whitelisted base
    out of the mask would leave a multi-base REF still overlapping it. The whole record's span
    is carved instead, which is what makes the whitelist mean "keep this variant"."""
    from plasgenomicsutils.lib import vcf_filters as F

    plain = tmp_path / "plain.bcf"
    F.tandem_repeat_mask(str(BCF), str(plain), bed=region_case["mask"])
    assert _positions(plain, region_case["picks"]) == []          # mask drops all four

    kept = tmp_path / "kept.bcf"
    n = F.tandem_repeat_mask(str(BCF), str(kept), bed=region_case["mask"],
                             keep_bed=region_case["keep"])
    assert n == 2
    assert _positions(kept, region_case["picks"]) == region_case["rescued"]


def test_a_whitelist_rescues_variants_outside_a_keep_region(tmp_path, region_case):
    """core_region_filter keeps what is inside its BED, so here the whitelist adds to it."""
    from plasgenomicsutils.lib import vcf_filters as F

    core = tmp_path / "core.bed"
    core.write_text(f"{region_case['chrom']}\t{max(region_case['picks']) + 10}"
                    f"\t{max(region_case['picks']) + 100000}\n")

    plain = tmp_path / "plain.bcf"
    F.core_region_filter(str(BCF), str(plain), bed=str(core))
    assert _positions(plain, region_case["picks"]) == []

    kept = tmp_path / "kept.bcf"
    n = F.core_region_filter(str(BCF), str(kept), bed=str(core), keep_bed=region_case["keep"])
    assert n == 2
    assert _positions(kept, region_case["picks"]) == region_case["rescued"]


def test_a_whitelist_does_not_keep_variants_the_region_rule_never_touched(tmp_path, region_case):
    """It exempts from this rule only -- it is not a general "keep these" instruction, and it
    must not resurrect anything the mask was not dropping in the first place."""
    from plasgenomicsutils.lib import vcf_filters as F

    with_wl, without = tmp_path / "a.bcf", tmp_path / "b.bcf"
    F.tandem_repeat_mask(str(BCF), str(with_wl), bed=region_case["mask"],
                         keep_bed=region_case["keep"])
    F.tandem_repeat_mask(str(BCF), str(without), bed=region_case["mask"])

    import subprocess

    def allpos(p):
        return set(subprocess.run(f"bcftools query -f '%POS\\n' {p}", shell=True,
                                  capture_output=True, text=True).stdout.split())
    # the only difference is the two rescued records
    assert allpos(with_wl) - allpos(without) == {str(p) for p in region_case["rescued"]}
    assert not allpos(without) - allpos(with_wl)


def test_a_whitelist_that_rescues_nothing_warns(tmp_path, region_case, capsys):
    """Silence here would be indistinguishable from a whitelist that worked, and the usual
    causes are a contig-name mismatch or 1-based positions written into a BED."""
    from plasgenomicsutils.lib import vcf_filters as F

    bad = tmp_path / "bad.bed"
    bad.write_text("nosuchchrom\t1\t2\n")
    n = F.tandem_repeat_mask(str(BCF), str(tmp_path / "o.bcf"), bed=region_case["mask"],
                             keep_bed=str(bad))
    assert n == 0
    assert "rescued nothing" in capsys.readouterr().out


def test_the_pipeline_threads_a_whitelist_through_a_region_step(tmp_path, region_case):
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    cfg = {"steps": [{"name": "tandem_repeat_mask",
                      "params": {"bed": region_case["mask"], "keep_bed": region_case["keep"]}}]}
    run_pipeline(str(BCF), str(tmp_path), cfg, emit_snp_bed=False)
    out = tmp_path / "01_tandem_repeat_mask.bcf"
    assert _positions(out, region_case["picks"]) == region_case["rescued"]


def test_one_bed_line_over_a_gene_rescues_every_variant_in_it(tmp_path):
    """The whitelist is region-based: a gene's interval keeps everything inside it, with no
    need to enumerate its variants."""
    import subprocess

    from plasgenomicsutils.lib import vcf_filters as F

    def positions(path, chrom):
        out = subprocess.run(f"bcftools query -f '%CHROM\\t%POS\\n' {path}", shell=True,
                             capture_output=True, text=True).stdout.strip().splitlines()
        return sorted(int(r.split("\t")[1]) for r in out if r.split("\t")[0] == chrom)

    chrom = subprocess.run(f"bcftools query -f '%CHROM\\n' {BCF} | head -1", shell=True,
                           capture_output=True, text=True).stdout.strip()
    start, end = 100_000, 130_000
    inside = [p for p in positions(str(BCF), chrom) if start < p <= end]
    assert len(inside) > 1                                # the region has to hold several

    mask = tmp_path / "mask.bed"
    mask.write_text(f"{chrom}\t0\t3000000\n")             # mask the whole chromosome
    gene = tmp_path / "gene.bed"
    gene.write_text(f"{chrom}\t{start}\t{end}\tsome_gene\n")   # one line, name column ignored

    bare, saved = tmp_path / "bare.bcf", tmp_path / "saved.bcf"
    F.tandem_repeat_mask(str(BCF), str(bare), bed=str(mask))
    n = F.tandem_repeat_mask(str(BCF), str(saved), bed=str(mask), keep_bed=str(gene))

    assert positions(bare, chrom) == []                   # all of it masked
    assert positions(saved, chrom) == inside              # exactly the gene's variants back
    assert n == len(inside)


# ---- the config is checked before anything runs -----------------------------------

def _cfg(**step):
    base = {"name": "maf_filter", "params": {"maf_min": 0.02}}
    base.update(step)
    return {"steps": [base]}


def test_a_pipeline_control_inside_params_is_named_as_such(tmp_path):
    """`enabled` in `params` reaches the step as an unexpected argument, which used to
    surface as a raw TypeError from inside the runner."""
    from plasgenomicsutils.lib.filter_pipeline import validate_config

    with pytest.raises(SystemExit, match='"enabled" goes on the step, not inside its params'):
        validate_config(_cfg(params={"maf_min": 0.02, "enabled": False}))
    # the same for the other keys that belong on the step
    for key in ("ext", "report", "name"):
        with pytest.raises(SystemExit, match=f'"{key}" goes on the step'):
            validate_config(_cfg(params={key: "x"}))


def test_an_unknown_param_names_what_the_step_takes(tmp_path):
    from plasgenomicsutils.lib.filter_pipeline import validate_config

    with pytest.raises(SystemExit) as e:
        validate_config(_cfg(params={"maf_mim": 0.02}))
    msg = str(e.value)
    assert "unknown param 'maf_mim'" in msg
    assert "Did you mean 'maf_min'?" in msg          # close enough to suggest
    assert "maf_filter takes:" in msg and "maf_max" in msg


def test_an_unknown_step_or_report_is_caught_here_too(tmp_path):
    from plasgenomicsutils.lib.filter_pipeline import validate_config

    with pytest.raises(SystemExit, match="Did you mean 'hard_qc_filter'"):
        validate_config(_cfg(name="hard_qc_filtr"))
    with pytest.raises(SystemExit, match="unknown pipeline report"):
        validate_config(_cfg(name="maf_filter", report=True))


def test_every_step_declares_what_it_accepts(tmp_path):
    """A step whose params cannot be introspected validates nothing, so the wrappers have to
    keep naming what they wrap."""
    from plasgenomicsutils.lib.filter_pipeline import (REPORTS, STEPS, _accepted_params)

    for name, fn in list(STEPS.items()) + list(REPORTS.items()):
        assert _accepted_params(fn) is not None, f"{name} accepts anything, so nothing is checked"


def test_the_default_config_passes_its_own_validation(tmp_path):
    from plasgenomicsutils.lib.filter_pipeline import DEFAULT_CONFIG, validate_config

    validate_config(DEFAULT_CONFIG)


def test_validation_happens_before_any_output_is_written(tmp_path):
    """The point of checking up front: a mistake in a late step must not cost the output of
    the early ones."""
    out = tmp_path / "run"
    cfg = {"steps": [{"name": "no_alt_filter"},
                     {"name": "maf_filter", "params": {"maf_mim": 0.02}}]}
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    with pytest.raises(SystemExit, match="unknown param"):
        run_pipeline(str(BCF), str(out), cfg, emit_snp_bed=False)
    assert not out.exists() or not list(out.iterdir())


def test_a_malformed_config_is_refused(tmp_path):
    from plasgenomicsutils.lib.filter_pipeline import validate_config

    with pytest.raises(SystemExit, match='non-empty "steps" list'):
        validate_config({"steps": []})
    with pytest.raises(SystemExit, match='has no "name"'):
        validate_config({"steps": [{"params": {}}]})
    with pytest.raises(SystemExit, match='"params" must be an object'):
        validate_config({"steps": [{"name": "maf_filter", "params": [1, 2]}]})


# ---- variant classes in the tally, and pruning intermediates -----------------------

def test_a_record_is_one_class_and_a_mixed_one_is_not_a_snp():
    from plasgenomicsutils.lib.bcftools import classify_record

    assert classify_record("A", "T") == "snps"
    assert classify_record("A", "T,G") == "snps"          # multiallelic, still all SNPs
    assert classify_record("A", "T,ATT") == "mixed"       # not counted as a SNP
    assert classify_record("ATT", "A") == "indels"
    assert classify_record("AT", "GC") == "mnps"
    assert classify_record("A", ".") == "no_alt"
    for symbolic in ("<*>", "<NON_REF>", "*", "A[chr1:100["):
        assert classify_record("A", symbolic) == "other"


def test_the_tally_carries_the_class_breakdown(tmp_path):
    from plasgenomicsutils.lib import filter_pipeline as P
    from plasgenomicsutils.lib.bcftools import VARIANT_TYPES

    cfg = {"steps": [{"name": "no_alt_filter"}]}
    tally = P.run_pipeline(str(BCF), str(tmp_path / "run"), cfg, emit_snp_bed=False)
    for row in tally:
        types = row["types"]
        # the classes account for every record, so the breakdown cannot drift from the count
        assert sum(types[t] for t in VARIANT_TYPES) == types["total"] == row["variants"]


def test_remove_intermediates_keeps_the_last_output_and_the_side_tables(tmp_path):
    from plasgenomicsutils.lib import filter_pipeline as P

    out = tmp_path / "run"
    cfg = {"remove_intermediates": True,
           "steps": [{"name": "no_alt_filter"}, {"name": "biallelic_snp_filter"},
                     {"name": "maf_filter", "params": {"maf_min": 0.0}}]}
    tally = P.run_pipeline(str(BCF), str(out), cfg, emit_snp_bed=False)
    left = sorted(p.name for p in out.iterdir())
    # the final callset, its index, and the record of how it was made -- nothing else
    assert left == ["03_maf_filter.bcf", "03_maf_filter.bcf.csi", "config_used.json"]
    # every step is still accounted for, and the removed ones say so
    steps = {r["step"]: r for r in tally}
    assert steps["no_alt_filter"]["removed"] is True
    assert "removed" not in steps["maf_filter"]
    assert steps["no_alt_filter"]["variants"] > 0        # counted before it was deleted


def test_the_input_is_never_removed(tmp_path):
    from plasgenomicsutils.lib import filter_pipeline as P

    """It is the caller's file, not an intermediate."""
    cfg = {"remove_intermediates": True, "steps": [{"name": "no_alt_filter"}]}
    P.run_pipeline(str(BCF), str(tmp_path / "run"), cfg, emit_snp_bed=False)
    assert BCF.exists()


def test_pruning_does_not_change_the_result(tmp_path):
    from plasgenomicsutils.lib import filter_pipeline as P

    """The only difference is what is left on disk."""
    cfg = {"steps": [{"name": "no_alt_filter"}, {"name": "biallelic_snp_filter"}]}
    a = P.run_pipeline(str(BCF), str(tmp_path / "keep"), cfg, emit_snp_bed=False)
    b = P.run_pipeline(str(BCF), str(tmp_path / "prune"), dict(cfg, remove_intermediates=True),
                       emit_snp_bed=False)
    import subprocess

    def positions(path):
        return subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\n", str(path)],
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                              text=True).stdout
    assert [r["variants"] for r in a] == [r["variants"] for r in b]
    assert positions(a[-1]["path"]) == positions(b[-1]["path"])


def test_an_unknown_top_level_key_is_caught(tmp_path):
    from plasgenomicsutils.lib.filter_pipeline import validate_config

    with pytest.raises(SystemExit, match="Did you mean 'remove_intermediates'"):
        validate_config({"remove_intermediate": True, "steps": [{"name": "maf_filter"}]})


# ---- the run records the cutoffs it used ------------------------------------------

def test_the_output_directory_records_the_config_as_run(tmp_path):
    """A config only holds what somebody wrote down; the thresholds that did the work live
    in the code. The run resolves them into its own directory."""
    import json

    from plasgenomicsutils.lib import filter_pipeline as P

    out = tmp_path / "run"
    P.run_pipeline(str(BCF), str(out), {"steps": [{"name": "maf_filter"}]},
                   emit_snp_bed=False)
    used = json.loads((out / "config_used.json").read_text())
    params = used["steps"][0]["params"]
    # every threshold the step actually applied, not just the one that was written down
    assert params["maf_min"] == 0.01 and params["maf_max"] is None
    assert set(params) >= {"maf_min", "maf_max", "meta", "group_col", "sample_col"}
    assert used["_meta"]["input"] == str(BCF.resolve())


def test_a_written_param_wins_over_the_default_it_records(tmp_path):
    import json

    from plasgenomicsutils.lib import filter_pipeline as P

    out = tmp_path / "run"
    P.run_pipeline(str(BCF), str(out),
                   {"steps": [{"name": "maf_filter", "params": {"maf_min": 0.05}}]},
                   emit_snp_bed=False)
    used = json.loads((out / "config_used.json").read_text())
    assert used["steps"][0]["params"]["maf_min"] == 0.05


def test_what_did_not_run_is_recorded_too(tmp_path):
    """A step switched off keeps its entry and its defaults -- the record is of the whole
    config, not only the parts that fired."""
    import json

    from plasgenomicsutils.lib import filter_pipeline as P

    out = tmp_path / "run"
    P.run_pipeline(str(BCF), str(out),
                   {"steps": [{"name": "no_alt_filter"},
                              {"name": "maf_filter", "enabled": False}]}, emit_snp_bed=False)
    used = json.loads((out / "config_used.json").read_text())
    off = used["steps"][1]
    assert off["enabled"] is False and off["params"]["maf_min"] == 0.01


def test_the_recorded_config_is_one_you_can_run_again(tmp_path):
    """It has to stay a valid config, or it is a log rather than a record."""
    import json

    from plasgenomicsutils.lib import filter_pipeline as P

    first = tmp_path / "a"
    P.run_pipeline(str(BCF), str(first), P.DEFAULT_CONFIG, emit_snp_bed=False)
    used = json.loads((first / "config_used.json").read_text())
    P.validate_config(used)                      # `_meta` and all
    again = tmp_path / "b"
    second = P.run_pipeline(str(BCF), str(again), used, emit_snp_bed=False)
    first_tally = json.loads((first / "config_used.json").read_text())
    assert first_tally["steps"] == used["steps"]
    assert [r.get("variants") for r in second] == [
        r.get("variants") for r in P.run_pipeline(str(BCF), str(tmp_path / "c"),
                                                  P.DEFAULT_CONFIG, emit_snp_bed=False)]


# ---- the MAF window means the same thing at any allele count -----------------------

def _multiallelic_vcf(tmp_path, n=50):
    """50 samples. 100: one ALT at 2%. 200: two ALTs at 2% each. 300: two at 1% each."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
           + "\t".join(f"s{i}" for i in range(n))]

    def row(pos, alt, carriers):
        gts, i = ["0/0"] * n, 0
        for idx, k in carriers.items():
            for _ in range(k):
                gts[i] = f"0/{idx}"
                i += 1
        return f"c1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts)

    hdr += [row(100, "T", {1: 2}), row(200, "T,G", {1: 2, 2: 2}),
            row(300, "T,G", {1: 1, 2: 1})]
    p = tmp_path / "multi.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _kept_positions(path):
    import subprocess
    out = subprocess.run(["bcftools", "query", "-f", "%POS\n", str(path)],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return [int(x) for x in out.stdout.split()]


def test_a_site_whose_every_allele_is_under_the_floor_is_dropped(tmp_path):
    """`bcftools -q` defaults to `nref`, the sum of the alternates: three at 1% each sum past
    a 2% floor. The floor is on the minor allele, so a site like that has to go."""
    from plasgenomicsutils.lib import vcf_filters as F

    out = str(tmp_path / "o.bcf")
    F.maf_filter(_multiallelic_vcf(tmp_path), out, maf_min=0.02)
    assert _kept_positions(out) == [100, 200]        # 300 is 0.01 + 0.01


def test_grouped_and_ungrouped_agree_on_a_multiallelic_site(tmp_path):
    """They compute the frequency by different routes -- `-q :minor` and `+fill-tags -t MAF`
    -- so a site they disagree on means one of the two is answering another question."""
    from plasgenomicsutils.lib import vcf_filters as F

    vcf = _multiallelic_vcf(tmp_path)
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tregion\n" + "".join(f"s{i}\tone\n" for i in range(50)))
    plain, grouped = str(tmp_path / "a.bcf"), str(tmp_path / "b.bcf")
    F.maf_filter(vcf, plain, maf_min=0.02)
    F.maf_filter(vcf, grouped, maf_min=0.02, meta=str(meta), group_col="region")
    assert _kept_positions(plain) == _kept_positions(grouped)


def test_a_symmetric_window_is_unchanged_by_the_switch(tmp_path):
    """On a biallelic callset the old `nref` window and `MAF >= floor` select the same
    records -- this fixed multiallelic sites rather than changing what the defaults mean."""
    import subprocess

    from plasgenomicsutils.lib import vcf_filters as F

    bi = str(tmp_path / "bi.bcf")
    subprocess.run(["bcftools", "view", "-m2", "-M2", "-v", "snps", str(BCF), "-Ob", "-o", bi],
                   check=True, stderr=subprocess.DEVNULL)
    tagged = str(tmp_path / "tagged.bcf")
    subprocess.run(["bcftools", "+fill-tags", bi, "-Ob", "-o", tagged, "--", "-t", "AC,AN,AF"],
                   check=True, stderr=subprocess.DEVNULL)
    old_window = subprocess.run(
        ["bcftools", "query", "-i", "AF>=0.02 && AF<=0.98", "-f", "%CHROM\t%POS\n", tagged],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout

    out = str(tmp_path / "o.bcf")
    F.maf_filter(bi, out, maf_min=0.02)
    got = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\n", out],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout
    assert old_window and got == old_window


def test_the_floor_is_the_second_most_common_allele_not_the_rarest(tmp_path):
    """`-q X:minor` would read the *least* frequent allele, which on 0.60/0.35/0.05 is 0.05 --
    dropping a site with a 35% minor allele. The floor has to be MAF."""
    from plasgenomicsutils.lib import vcf_filters as F

    n = 100
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
           + "\t".join(f"s{i}" for i in range(n))]
    gts = ["1/1"] * 35 + ["2/2"] * 5 + ["0/0"] * 60       # REF .60, ALT1 .35, ALT2 .05
    hdr.append("c1\t100\t.\tA\tT,G\t.\t.\t.\tGT\t" + "\t".join(gts))
    v = tmp_path / "tri.vcf"
    v.write_text("\n".join(hdr) + "\n")

    out = str(tmp_path / "o.bcf")
    F.maf_filter(str(v), out, maf_min=0.10)
    assert _kept_positions(out) == [100]


def test_the_expression_reads_maf_and_only_bounds_the_top_when_asked():
    from plasgenomicsutils.lib.vcf_filters import _maf_expr

    assert _maf_expr(0.02, None) == "MAF >= 0.02"
    # an asymmetric window is a band on the alternate, not a minor-allele floor: otherwise a
    # site at AF 0.75 would fail a 0.3 floor while sitting inside a 0.3-0.8 band
    assert _maf_expr(0.3, 0.8) == "MAX(AF) >= 0.3 && MAX(AF) <= 0.8"
