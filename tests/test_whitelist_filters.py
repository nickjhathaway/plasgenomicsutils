"""Whitelisting (`keep_bed`) on the filters that judge a record on its own numbers.

The region filters fold a whitelist into their BED and are covered in
`test_vcf_region_and_biallelic.py`; these cover the ones that cannot -- a QC metric, an
allele count, a missingness rate, an allele frequency.

They shell out to bcftools/bedtools and skip cleanly where those aren't installed.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import filter_pipeline as P
from plasgenomicsutils.lib import vcf_filters as F
from plasgenomicsutils.lib.bcftools import q, sh

needs_both = pytest.mark.skipif(
    not (shutil.which("bcftools") and shutil.which("bedtools")),
    reason="bcftools+bedtools not on PATH")

pytestmark = needs_both

def _wl_vcf(tmp_path, n_samples=40):
    """3 records: common/clean, a rare low-QD one, and a rare-or-low-QD control."""
    samples = [f"s{i}" for i in range(1, n_samples + 1)]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=Pf3D7_07_v3,length=1445207>"]
    for t in ("QD", "MQ", "SOR", "MQRankSum", "ReadPosRankSum"):
        hdr.append(f'##INFO=<ID={t},Number=1,Type=Float,Description="{t}">')
    hdr += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]

    def rec(pos, n_alt, qd, name):
        gts = ["1/1"] * n_alt + ["0/0"] * (n_samples - n_alt)
        fmt = "\t".join(f"{g}:{'0,30' if g == '1/1' else '30,0'}" for g in gts)
        info = f"QD={qd};MQ=60;SOR=1.0;MQRankSum=0.0;ReadPosRankSum=0.0"
        return f"Pf3D7_07_v3\t{pos}\t{name}\tA\tT\t100\t.\t{info}\tGT:AD\t{fmt}"

    vcf = tmp_path / "in.vcf"
    vcf.write_text("\n".join(hdr + [rec(1000, 20, 30, "clean"),
                                    rec(2000, 1, 5, "resistance"),
                                    rec(3000, 1, 5, "control")]) + "\n")
    bcf = tmp_path / "in.bcf"
    sh(f"bcftools view {q(str(vcf))} -Ob -o {q(str(bcf))}")
    sh(f"bcftools index -f {q(str(bcf))}")
    keep = tmp_path / "keep.bed"
    keep.write_text("Pf3D7_07_v3\t1999\t2000\tresistance\n")   # 0-based half-open
    return str(bcf), str(keep)


def _ids(path):
    p = subprocess.run(["bcftools", "query", "-f", "%ID\n", str(path)],
                       stdout=subprocess.PIPE, text=True)
    return p.stdout.split()


@pytest.mark.parametrize("fn,kwargs", [
    (F.hard_qc_filter, {}),
    (F.singleton_add_ads, {}),
    (F.maf_filter, {"maf_min": 0.05}),
])
def test_keep_bed_rescues_only_what_it_covers(tmp_path, fn, kwargs):
    inp, keep = _wl_vcf(tmp_path)
    plain, wl = str(tmp_path / "plain.bcf"), str(tmp_path / "wl.bcf")

    fn(inp, plain, **kwargs)
    assert _ids(plain) == ["clean"]                       # both rare/failing ones dropped

    rescued = fn(inp, wl, keep_bed=keep, **kwargs)
    assert rescued == 1
    assert _ids(wl) == ["clean", "resistance"]            # the control is still gone


def test_keep_bed_leaves_the_unwhitelisted_run_untouched(tmp_path):
    """Passing no whitelist must go down the original single-pipe path unchanged."""
    inp, _ = _wl_vcf(tmp_path)
    a, b = str(tmp_path / "a.bcf"), str(tmp_path / "b.bcf")
    F.maf_filter(inp, a, maf_min=0.05)
    F.maf_filter(inp, b, maf_min=0.05, keep_bed=None)
    assert _ids(a) == _ids(b) == ["clean"]


def test_rescued_singleton_record_carries_the_added_ads_tag(tmp_path):
    """The rescue reads from the tagged intermediate, or concat would reject the header."""
    inp, keep = _wl_vcf(tmp_path)
    out = str(tmp_path / "out.bcf")
    F.singleton_add_ads(inp, out, keep_bed=keep)
    p = subprocess.run(["bcftools", "query", "-f", "%ID\t[%ADS ]\n", out],
                       stdout=subprocess.PIPE, text=True)
    rows = dict(line.split("\t", 1) for line in p.stdout.strip().splitlines())
    assert set(rows) == {"clean", "resistance"}
    assert all(vals.strip() and "." not in vals.split() for vals in rows.values())


def test_pipeline_keep_bed_applies_to_every_whitelistable_step(tmp_path):
    inp, keep = _wl_vcf(tmp_path)
    cfg = {"keep_bed": keep,
           "steps": [{"name": "hard_qc_filter"},
                     {"name": "singleton_filter_add_ads"},
                     {"name": "maf_filter", "params": {"maf_min": 0.05}}]}
    tally = P.run_pipeline(inp, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    steps = {r["step"]: r for r in tally if r["step"] != "input"}
    assert all(steps[s].get("rescued") == 1 for s in
               ("hard_qc_filter", "singleton_filter_add_ads", "maf_filter"))
    assert "resistance" in _ids(steps["maf_filter"]["path"])

    # without it the resistance variant is lost at the first step it fails
    cfg.pop("keep_bed")
    tally = P.run_pipeline(inp, str(tmp_path / "run2"), cfg, emit_snp_bed=False)
    last = [r for r in tally if r["step"] == "maf_filter"][0]
    assert "resistance" not in _ids(last["path"])
    assert not last.get("rescued")


def test_a_step_opts_out_of_the_pipeline_wide_whitelist(tmp_path):
    inp, keep = _wl_vcf(tmp_path)
    cfg = {"keep_bed": keep,
           "steps": [{"name": "maf_filter",
                      "params": {"maf_min": 0.05, "keep_bed": None}}]}
    tally = P.run_pipeline(inp, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    last = [r for r in tally if r["step"] == "maf_filter"][0]
    assert "resistance" not in _ids(last["path"])


def test_a_whitelist_that_works_somewhere_does_not_warn_anywhere(tmp_path, capsys):
    """The warning is a verdict on the run, not on each step: over a chain, most steps have
    nothing for a whitelist to do, and saying so at every one reads as a misconfiguration."""
    inp, keep = _wl_vcf(tmp_path)
    cfg = {"keep_bed": keep,
           "steps": [{"name": "no_alt_filter"},          # nothing here to rescue
                     {"name": "hard_qc_filter"}]}        # rescues the resistance variant
    tally = P.run_pipeline(inp, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    steps = {r["step"]: r for r in tally if r["step"] != "input"}
    assert not steps["no_alt_filter"].get("rescued")
    assert steps["hard_qc_filter"]["rescued"] == 1
    assert "rescued nothing" not in capsys.readouterr().out


def test_a_whitelist_that_matches_nothing_warns_once_for_the_whole_run(tmp_path, capsys):
    inp, _ = _wl_vcf(tmp_path)
    bad = tmp_path / "bad.bed"
    bad.write_text("nosuchchrom\t1\t2\n")
    cfg = {"keep_bed": str(bad),
           "steps": [{"name": "no_alt_filter"},
                     {"name": "hard_qc_filter"},
                     {"name": "maf_filter", "params": {"maf_min": 0.05}}]}
    P.run_pipeline(inp, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    out = capsys.readouterr().out
    assert out.count("rescued nothing") == 1
    assert "3 step(s)" in out


def test_a_single_filter_run_on_its_own_still_warns_immediately(tmp_path, capsys):
    """Outside a pipeline there is no later step to redeem it, so the warning is the same
    verdict, reported where it happens."""
    inp, _ = _wl_vcf(tmp_path)
    bad = tmp_path / "bad.bed"
    bad.write_text("nosuchchrom\t1\t2\n")
    F.hard_qc_filter(inp, str(tmp_path / "o.bcf"), keep_bed=str(bad))
    assert "rescued nothing" in capsys.readouterr().out


def test_the_biallelic_filter_names_the_whitelisted_records_it_removes(tmp_path, capsys):
    """A whitelist that silently does nothing reads as a guarantee it is not making.

    `biallelic_snp_filter` deliberately does not honour `keep_bed` (decision 3): letting one
    multiallelic record past `-M2` because it sits under a whitelist would break anything
    downstream relying on the file being biallelic, silently and far from here. So the
    whitelist becomes a diagnostic -- it says what it would have saved.
    """
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    hdr = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n")
    body = ("chr1\t1000\t.\tA\tC,G\t.\t.\t.\tGT:AD\t1/1:0,9,0\t1/1:0,9,0\t2/2:0,0,9\t"
            "2/2:0,0,9\n"
            "chr1\t2000\t.\tA\tC\t.\t.\t.\tGT:AD\t1/1:0,9\t1/1:0,9\t0/0:9,0\t0/0:9,0\n")
    src = tmp_path / "in.vcf"
    src.write_text(hdr + body)
    bed = tmp_path / "wl.bed"
    bed.write_text("chr1\t999\t1001\n")          # covers the multiallelic record only

    cfg = {"keep_bed": str(bed),
           "steps": [{"name": "biallelic_snp_filter",
                      "params": {"snps_only": True, "biallelic": True,
                                 "mnp_handling": "remove"}}]}
    run_pipeline(str(src), str(tmp_path / "out"), cfg, emit_snp_bed=False)
    printed = capsys.readouterr().out
    assert "whitelisted record(s) were removed by biallelic_snp_filter" in printed
    assert "chr1:999 A>C,G" in printed


def test_nothing_is_said_when_the_whitelist_loses_nothing(tmp_path, capsys):
    from plasgenomicsutils.lib.filter_pipeline import run_pipeline

    hdr = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n")
    body = "chr1\t2000\t.\tA\tC\t.\t.\t.\tGT:AD\t1/1:0,9\t1/1:0,9\t0/0:9,0\t0/0:9,0\n"
    src = tmp_path / "in.vcf"
    src.write_text(hdr + body)
    bed = tmp_path / "wl.bed"
    bed.write_text("chr1\t1999\t2001\n")
    cfg = {"keep_bed": str(bed),
           "steps": [{"name": "biallelic_snp_filter",
                      "params": {"snps_only": True, "biallelic": True,
                                 "mnp_handling": "remove"}}]}
    run_pipeline(str(src), str(tmp_path / "out"), cfg, emit_snp_bed=False)
    assert "were removed by biallelic_snp_filter" not in capsys.readouterr().out
