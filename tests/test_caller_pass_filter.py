"""The caller's FILTER column is its own question, asked in its own step.

A joint callset arrives with the caller's verdicts already in FILTER -- Pf7's VQSR tranches,
`Low_VQSLOD`, its region classes, `MissingVQSLOD` on the contigs VQSR never scored. Until
v0.3.2 `hard_qc_filter` enforced those silently, by selecting PASS at the end: a
mitochondrial subset whose every record was `MissingVQSLOD;Mitochondrion` came out empty
with no word as to why. Now `caller_pass_filter` acts on the FILTER column and says what it
removed, and `hard_qc_filter` judges on its metrics alone.
"""

import pathlib
import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import vcf_filters as F
from plasgenomicsutils.lib import filter_pipeline as P

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not on PATH")

DATA = pathlib.Path(__file__).parent / "data"
BCF = DATA / "ghana_cambodia.pf7.tiny.bcf"

GOOD = "QD=30;MQ=60;SOR=1;MQRankSum=0;ReadPosRankSum=0"
BAD = "QD=5;MQ=60;SOR=1;MQRankSum=0;ReadPosRankSum=0"


def _vcf(tmp_path, rows, name="in.vcf"):
    """(pos, FILTER, INFO) rows in a GATK-shaped callset with the FILTER ids declared."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>"]
    for t in ("QD", "MQ", "SOR", "MQRankSum", "ReadPosRankSum"):
        hdr.append(f'##INFO=<ID={t},Number=1,Type=Float,Description="{t}">')
    for fid in ("Low_VQSLOD", "MissingVQSLOD", "Mitochondrion", "SubtelomericRepeat"):
        hdr.append(f'##FILTER=<ID={fid},Description="{fid}">')
    hdr.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
    hdr.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1")
    for pos, flt, info in rows:
        hdr.append(f"chr1\t{pos}\t.\tA\tT\t222\t{flt}\t{info}\tGT\t0/1")
    p = tmp_path / name
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _kept(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\t%FILTER\n", path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return {int(a): b for a, b in (ln.split("\t") for ln in out.stdout.splitlines())}


ROWS = [(1000, ".", GOOD), (2000, "PASS", GOOD), (3000, "Low_VQSLOD", GOOD),
        (4000, "MissingVQSLOD;Mitochondrion", GOOD), (5000, "PASS", BAD),
        (6000, "Low_VQSLOD", BAD)]


# ---- hard_qc_filter no longer enforces the caller's flags ----------------------------

def test_hard_qc_judges_on_its_metrics_and_leaves_the_callers_flags_alone(tmp_path, capsys):
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(_vcf(tmp_path, ROWS), out)
    kept = _kept(out)
    assert sorted(kept) == [1000, 2000, 3000, 4000]         # the two QD=5 records are gone
    # the caller's flags survive on the kept records, and the step says so
    assert kept[3000] == "Low_VQSLOD" and kept[4000] == "MissingVQSLOD;Mitochondrion"
    assert kept[1000] == "PASS"                              # bcftools' own PASS on a clean '.'
    cap = capsys.readouterr()
    err = cap.out + cap.err
    assert "2 kept record(s) carry a FILTER set by the caller" in err
    assert "Low_VQSLOD: 1" in err and "caller_pass_filter" in err


def test_a_callset_the_caller_flagged_entirely_is_not_emptied_by_hard_qc(tmp_path):
    """The case that prompted this: a mitochondrial subset, every record
    MissingVQSLOD;Mitochondrion, 6 of 19 clean on the metrics -- and 0 came out."""
    rows = [(p, "MissingVQSLOD;Mitochondrion", GOOD if p % 2 else BAD) for p in range(1, 20)]
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(_vcf(tmp_path, rows), out)
    assert len(_kept(out)) == 10


def test_a_whitelisted_record_still_carries_fail_and_nothing_else_is_said_about_it(tmp_path, capsys):
    bed = tmp_path / "keep.bed"
    bed.write_text("chr1\t4999\t5000\n")
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(_vcf(tmp_path, ROWS), out, keep_bed=str(bed))
    kept = _kept(out)
    assert kept[5000] == "FAIL"                              # rescued, and says it failed
    # FAIL is this package's verdict, not the caller's, so it is not in the caller-flag note
    cap = capsys.readouterr()
    assert "caller (Low_VQSLOD: 1, MissingVQSLOD;Mitochondrion: 1)" in cap.out + cap.err


# ---- caller_pass_filter ------------------------------------------------------------------

def test_pass_and_missing_are_kept_and_every_flag_removes(tmp_path, capsys):
    out = str(tmp_path / "o.bcf")
    F.caller_pass_filter(_vcf(tmp_path, ROWS), out)
    assert sorted(_kept(out)) == [1000, 2000, 5000]          # the metrics are not its business
    cap = capsys.readouterr()
    err = cap.out + cap.err
    assert "3 record(s) removed on the caller's own FILTER" in err
    assert "Low_VQSLOD: 2" in err and "MissingVQSLOD;Mitochondrion: 1" in err


def test_allow_tolerates_flags_but_only_when_every_flag_on_the_record_is_allowed(tmp_path):
    a, b = str(tmp_path / "a.bcf"), str(tmp_path / "b.bcf")
    F.caller_pass_filter(_vcf(tmp_path, ROWS), a, allow=["MissingVQSLOD", "Mitochondrion"])
    assert sorted(_kept(a)) == [1000, 2000, 4000, 5000]
    # allowing one of the two flags on 4000 is not enough
    F.caller_pass_filter(_vcf(tmp_path, ROWS, name="in2.vcf"), b, allow=["MissingVQSLOD"])
    assert sorted(_kept(b)) == [1000, 2000, 5000]


def test_a_header_with_nothing_to_act_on_copies_through_and_says_so(tmp_path, capsys):
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
           "chr1\t1000\t.\tA\tT\t222\t.\t.\tGT\t0/1"]
    v = tmp_path / "bare.vcf"
    v.write_text("\n".join(hdr) + "\n")
    out = str(tmp_path / "o.bcf")
    F.caller_pass_filter(str(v), out)
    assert sorted(_kept(out)) == [1000]
    cap = capsys.readouterr()
    assert "declares no FILTER this step acts on" in cap.out + cap.err


def test_the_whitelist_rescues_a_flagged_record_with_its_flag(tmp_path):
    bed = tmp_path / "keep.bed"
    bed.write_text("chr1\t2999\t3000\n")
    out = str(tmp_path / "o.bcf")
    assert F.caller_pass_filter(_vcf(tmp_path, ROWS), out, keep_bed=str(bed)) == 1
    assert _kept(out)[3000] == "Low_VQSLOD"


def test_the_default_chain_runs_it_before_hard_qc_and_it_is_whitelistable():
    names = [s["name"] for s in P.DEFAULT_CONFIG["steps"]]
    assert names.index("caller_pass_filter") < names.index("hard_qc_filter")
    assert "caller_pass_filter" in P.WHITELISTABLE and "caller_pass_filter" in P.STEPS


def test_on_the_pf7_fixture_the_split_reproduces_what_pass_selection_used_to_do(tmp_path):
    """Before: hard_qc_filter kept 'PASS after the metrics'. Now the same set is the two
    steps in sequence, and the FILTER-only removals are their own number."""
    a = str(tmp_path / "a.bcf")
    F.caller_pass_filter(str(BCF), a)
    b = str(tmp_path / "b.bcf")
    F.hard_qc_filter(a, b)
    both = len(_kept(b))
    # the old single step, reproduced by hand
    ref = str(tmp_path / "ref.bcf")
    subprocess.run("bcftools filter -m + -s FAIL -e 'QD < 10 || MQ < 55 || SOR > 3 || "
                   "(MQRankSum!=\".\" && MQRankSum < -5) || "
                   "(ReadPosRankSum!=\".\" && ReadPosRankSum < -5)' "
                   f"{BCF} -Ou | bcftools view -f PASS -Ob -o {ref}", shell=True, check=True)
    assert both == len(_kept(ref)) > 0


# ---- an emptied callset stops the chain cleanly ------------------------------------------

def test_a_step_that_removes_everything_stops_the_run_and_records_the_rest_as_not_run(tmp_path, capsys):
    cfg = {"steps": [
        {"name": "no_alt_filter", "params": {"trim": True}},
        {"name": "hard_qc_filter", "params": {"qd": 1e9}},       # nothing can pass this
        {"name": "singleton_counts", "report": True, "ext": "tsv"},
        {"name": "singleton_filter_add_ads"},
        {"name": "sample_coverage_filter"},
        {"name": "locus_missingness_filter"},
    ]}
    tally = P.run_pipeline(str(BCF), str(tmp_path), cfg)      # no exception
    steps = {t["step"]: t for t in tally[1:]}
    assert steps["hard_qc_filter"]["variants"] == 0
    for name in ("singleton_counts", "singleton_filter_add_ads", "sample_coverage_filter",
                 "locus_missingness_filter"):
        assert steps[name] == {"step": name, "skipped": True, "reason": "no variants remain"}
    assert "snp_bed" not in steps                              # no panel of nothing
    cap = capsys.readouterr()
    assert "no variants remain after [02] hard_qc_filter; 4 later step(s) not run" in cap.out + cap.err
    assert not (tmp_path / "04_singleton_filter_add_ads.bcf").exists()


def test_sample_coverage_over_no_loci_drops_nobody(tmp_path, capsys):
    empty = str(tmp_path / "empty.bcf")
    F.singleton_add_ads(str(BCF), str(tmp_path / "ads.bcf"))
    subprocess.run(f"bcftools view -i 'QUAL<0' {tmp_path}/ads.bcf -Ob -o {empty}",
                   shell=True, check=True)
    out = str(tmp_path / "o.bcf")
    assert F.sample_coverage_filter(empty, out) == []
    cap = capsys.readouterr()
    assert "no variants to measure coverage on; no sample dropped" in cap.out + cap.err
    n = subprocess.run(["bcftools", "query", "-l", out], capture_output=True, text=True)
    assert len(n.stdout.split()) == 60                         # every sample still there


def test_locus_missingness_without_ads_is_refused_by_name(tmp_path):
    with pytest.raises(SystemExit, match="no FORMAT/ADS.*singleton_filter_add_ads"):
        F.locus_missingness_filter(str(BCF), str(tmp_path / "o.bcf"))


# ---- reset_filter: the caller's column cleared before the chain runs ----------------------

def test_reset_filter_clears_every_flag_and_says_what_it_cleared(tmp_path, capsys):
    out = str(tmp_path / "o.bcf")
    cleared = F.reset_filter(_vcf(tmp_path, ROWS), out)
    assert cleared == {"Low_VQSLOD": 2, "MissingVQSLOD;Mitochondrion": 1}
    assert set(_kept(out).values()) == {"."}
    cap = capsys.readouterr()
    assert "3 carried a caller flag (Low_VQSLOD: 2, MissingVQSLOD;Mitochondrion: 1)" in cap.out + cap.err


def test_the_config_key_runs_it_as_step_zero_and_the_chain_then_sees_no_flags(tmp_path, capsys):
    v = _vcf(tmp_path, ROWS)
    cfg = {"reset_filter": True, "steps": [{"name": "caller_pass_filter"},
                                          {"name": "hard_qc_filter"}]}
    tally = P.run_pipeline(v, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    assert [t["step"] for t in tally] == ["input", "reset_filter", "caller_pass_filter",
                                          "hard_qc_filter"]
    assert tally[1]["variants"] == 6 and tally[1]["path"].endswith("00_reset_filter.bcf")
    assert tally[2]["variants"] == 6                          # nothing left for it to act on
    assert tally[3]["variants"] == 4                          # the metrics alone decide
    cap = capsys.readouterr()
    assert "[00] reset_filter" in cap.out + cap.err
    assert "declares no FILTER this step acts on" in cap.out + cap.err
    import json
    used = json.loads((tmp_path / "run" / "config_used.json").read_text())
    assert used["reset_filter"] is True


def test_reset_filter_is_off_by_default_and_recorded_as_such(tmp_path):
    import json
    assert P.DEFAULT_CONFIG["reset_filter"] is False
    P.validate_config(P.DEFAULT_CONFIG)
    v = _vcf(tmp_path, ROWS)
    P.run_pipeline(v, str(tmp_path / "run"), {"steps": [{"name": "caller_pass_filter"}]},
                   emit_snp_bed=False)
    used = json.loads((tmp_path / "run" / "config_used.json").read_text())
    assert used["reset_filter"] is False
    assert not (tmp_path / "run" / "00_reset_filter.bcf").exists()


def test_the_cli_flag_sets_it(tmp_path, monkeypatch):
    import json
    from plasgenomicsutils.scripts.vcf import filter_pipeline as CLI
    v = _vcf(tmp_path, ROWS)
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"steps": [{"name": "caller_pass_filter"}]}))
    monkeypatch.setattr("sys.argv", ["filter_pipeline", "--input", v, "--config", str(cfg),
                                     "--outdir", str(tmp_path / "run"), "--reset-filter",
                                     "--no-snp-bed"])
    CLI.filter_pipeline()
    used = json.loads((tmp_path / "run" / "config_used.json").read_text())
    assert used["reset_filter"] is True
    assert (tmp_path / "run" / "00_reset_filter.bcf").exists()
    counts = (tmp_path / "run" / "variant_counts.tsv").read_text()
    assert "reset_filter" in counts
