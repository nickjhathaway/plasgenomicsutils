"""The filter pipeline's report steps: diagnostics that read the callset but never
change it, and the ordering trap that makes one of them meaningless."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"
BCF = DATA / "ghana_cambodia.pf7.tiny.bcf"


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
    # every other step runs
    assert all(s.get("enabled", True) for s in DEFAULT_CONFIG["steps"]
               if s["name"] != "paralog_mask")


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
    assert rows[0] == "step\tkind\tcount\tpath"
    kinds = {r.split("\t")[0]: r.split("\t")[1] for r in rows[1:]}
    assert kinds["singleton_counts"] == "report"
    assert kinds["paralog_mask"] == "skipped"
    assert kinds["singleton_filter_add_ads"] == "variants"
