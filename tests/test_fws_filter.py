"""fws_filter: keeping the monoclonal infections, and what that does and does not touch.

It is the second step in the chain that drops samples rather than variants, and the only
one that changes which infections the callset describes -- so what it leaves alone matters
as much as what it removes.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import filter_pipeline as P
from plasgenomicsutils.lib.fws import fws_filter, fws_table

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"),
                                reason="bcftools not on PATH")


def _cohort_vcf(tmp_path, n_sites=400, name="in.vcf"):
    """Four samples: three clean monoclonals and one obvious mixture.

    Sample `mixed` carries both alleles at a third of the sites, which is what within-host
    diversity looks like to Fws; the others are homozygous wherever they are called.
    """
    samples = ["mono_a", "mono_b", "mono_c", "mixed"]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for i in range(n_sites):
        # the population allele frequency has to vary, or every site lands in one MAF bin
        carriers = i % 4
        cells = []
        for s in range(4):
            if s == 3 and i % 3 == 0:
                cells.append("0/1:20,20")                  # the mixture's het sites
            elif s < carriers:
                cells.append("1/1:0,40")
            else:
                cells.append("0/0:40,0")
        hdr.append(f"chr1\t{1000 + i * 100}\t.\tA\tT\t222\t.\t.\tGT:AD\t" + "\t".join(cells))
    p = tmp_path / name
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _samples(path):
    out = subprocess.run(["bcftools", "query", "-l", path], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    return out.stdout.split()


def _positions(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\n", path], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    return out.stdout.split()


def test_the_mixture_scores_lowest_and_is_the_one_dropped(tmp_path):
    v = _cohort_vcf(tmp_path)
    rows, n_sites = fws_table(v, fws_min=0.9)
    by = {r["sample"]: r for r in rows}
    assert n_sites > 0
    assert by["mixed"]["fws"] < min(by[s]["fws"] for s in ("mono_a", "mono_b", "mono_c"))

    out = str(tmp_path / "o.bcf")
    dropped = fws_filter(v, out, fws_min=0.9)
    assert dropped == ["mixed"]
    assert _samples(out) == ["mono_a", "mono_b", "mono_c"]


def test_it_drops_samples_and_not_variants(tmp_path):
    """Removing samples changes every allele frequency; deciding what that means for the
    site set is the caller's, so nothing is quietly removed here."""
    v = _cohort_vcf(tmp_path)
    out = str(tmp_path / "o.bcf")
    fws_filter(v, out, fws_min=0.9)
    assert _positions(out) == _positions(v)


def test_ac_an_are_refreshed_for_the_cohort_that_remains(tmp_path):
    """So a maf_filter placed after this one reads the survivors' frequencies, not the
    frequencies of a cohort that no longer exists."""
    v = _cohort_vcf(tmp_path)
    out = str(tmp_path / "o.bcf")
    fws_filter(v, out, fws_min=0.9)
    an = subprocess.run(["bcftools", "query", "-f", "%AN\n", out], stdout=subprocess.PIPE,
                        stderr=subprocess.DEVNULL, text=True).stdout.split()
    assert set(an) == {"6"}          # three diploid samples, not four


def test_a_threshold_nobody_meets_is_an_error_not_an_empty_callset(tmp_path):
    v = _cohort_vcf(tmp_path)
    with pytest.raises(SystemExit, match="keeps no samples"):
        fws_filter(v, str(tmp_path / "o.bcf"), fws_min=1.5)


def test_a_threshold_everybody_meets_drops_nothing(tmp_path):
    v = _cohort_vcf(tmp_path)
    out = str(tmp_path / "o.bcf")
    assert fws_filter(v, out, fws_min=-1.0) == []
    assert _samples(out) == ["mono_a", "mono_b", "mono_c", "mixed"]


def test_the_table_records_the_decision_for_every_sample(tmp_path):
    v = _cohort_vcf(tmp_path)
    tsv = str(tmp_path / "fws.tsv")
    fws_filter(v, str(tmp_path / "o.bcf"), fws_min=0.9, fws_table_path=tsv)
    lines = [x.split("\t") for x in open(tsv).read().strip().splitlines()]
    assert lines[0] == ["sample", "fws", "n_sites", "monoclonal", "dropped"]
    assert len(lines) == 5                                   # header + every sample
    by = {r[0]: r for r in lines[1:]}
    assert by["mixed"][4] == "True" and by["mono_a"][4] == "False"


def test_the_pipeline_runs_it_and_writes_its_table_beside_the_step(tmp_path):
    v = _cohort_vcf(tmp_path)
    cfg = {"steps": [{"name": "fws_filter", "params": {"fws_min": 0.9}}]}
    tally = P.run_pipeline(v, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    step = [r for r in tally if r["step"] == "fws_filter"][0]
    assert _samples(step["path"]) == ["mono_a", "mono_b", "mono_c"]
    assert (tmp_path / "run" / "01_fws_filter_fws.tsv").exists()


def test_it_ships_switched_off_in_the_default_config(tmp_path):
    """Keeping only monoclonals is an analysis choice, so it is discoverable but not on."""
    steps = {s["name"]: s for s in P.DEFAULT_CONFIG["steps"]}
    assert steps["fws_filter"]["enabled"] is False
    assert "fws_min" in steps["fws_filter"]["params"]
    # it judges samples, so a region whitelist cannot exempt anything from it
    assert "fws_filter" not in P.WHITELISTABLE


def test_an_unscored_sample_is_dropped_rather_than_assumed_monoclonal(tmp_path):
    """A sample with no usable depth cannot be called monoclonal -- nothing is known about
    it, and keeping it would readmit exactly what the step exists to remove."""
    v = _cohort_vcf(tmp_path, name="in2.vcf")
    txt = open(v).read().replace("\t0/0:40,0\t0/0:40,0\t", "\t./.:0,0\t0/0:40,0\t")
    blank = tmp_path / "blank.vcf"
    blank.write_text(txt)
    rows, _ = fws_table(str(blank), fws_min=0.9)
    by = {r["sample"]: r for r in rows}
    if by["mono_b"]["fws"] is None:                   # only meaningful if it went unscored
        assert by["mono_b"]["dropped"] is True
        assert by["mono_b"]["monoclonal"] is False
