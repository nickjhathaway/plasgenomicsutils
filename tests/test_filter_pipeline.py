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
