"""Cohort statistics judged on a named subset, with every sample's genotypes kept.

A frequency, a missingness rate and a carrier count are properties of *who is in the file*,
so the same locus passes in an analysis cohort and fails in a larger callset that merely
dilutes it. `stat_samples` lets a superset callset be filtered as though it were the cohort,
which is what stops a 249-sample analysis losing loci when it is run inside a 374-sample
call. Nothing is removed from the output but records: every sample's genotypes survive.
"""

import pathlib
import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import filter_pipeline as P
from plasgenomicsutils.lib import vcf_filters as F

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not on PATH")


def _vcf(tmp_path, rows, n_cohort=10, n_extra=30, name="in.vcf"):
    """`rows` are (pos, n_alt_cohort, n_alt_extra): how many samples of each group carry
    the alternate. Cohort samples are c00.., extras e00..."""
    cohort = [f"c{i:02d}" for i in range(n_cohort)]
    extra = [f"e{i:02d}" for i in range(n_extra)]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(cohort + extra)]
    for pos, n_c, n_e in rows:
        cols = []
        for i in range(n_cohort):
            cols.append("1/1:0,30" if i < n_c else "0/0:30,0")
        for i in range(n_extra):
            cols.append("1/1:0,30" if i < n_e else "0/0:30,0")
        hdr.append(f"chr1\t{pos}\t.\tA\tT\t500\t.\t.\tGT:AD\t" + "\t".join(cols))
    p = tmp_path / name
    p.write_text("\n".join(hdr) + "\n")
    return str(p), cohort, extra


def _kept(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\n", path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return [int(x) for x in out.stdout.split()]


def _samples(path):
    out = subprocess.run(["bcftools", "query", "-l", path], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    return out.stdout.split()


def _list_file(tmp_path, names, name="cohort.txt"):
    p = tmp_path / name
    p.write_text("\n".join(names) + "\n")
    return str(p)


# ---- maf_filter ---------------------------------------------------------------------

def test_a_frequency_floor_can_be_judged_on_the_cohort_inside_a_bigger_callset(tmp_path):
    """Position 1000 is at 30% in the cohort and 7.5% overall; 2000 is the reverse. A 20%
    floor keeps the first and drops the second when the cohort decides."""
    v, cohort, extra = _vcf(tmp_path, [(1000, 3, 0), (2000, 0, 12)])
    plain = str(tmp_path / "plain.bcf")
    F.maf_filter(v, plain, maf_min=0.20)
    assert _kept(plain) == [2000]                       # the whole callset's frequencies
    sub = str(tmp_path / "sub.bcf")
    F.maf_filter(v, sub, maf_min=0.20, stat_samples=_list_file(tmp_path, cohort))
    assert _kept(sub) == [1000]
    assert _samples(sub) == cohort + extra              # every sample's genotypes kept


def test_the_genotypes_are_untouched_only_the_verdict_changes(tmp_path):
    v, cohort, _ = _vcf(tmp_path, [(1000, 3, 0)])
    sub = str(tmp_path / "sub.bcf")
    F.maf_filter(v, sub, maf_min=0.20, stat_samples=_list_file(tmp_path, cohort))
    gt = subprocess.run(["bcftools", "query", "-f", "[%GT ]\n", sub], stdout=subprocess.PIPE,
                        stderr=subprocess.DEVNULL, text=True).stdout.split()
    assert gt.count("1/1") == 3 and gt.count("0/0") == 37


def test_a_comma_separated_list_works_as_well_as_a_file(tmp_path):
    v, cohort, _ = _vcf(tmp_path, [(1000, 3, 0)])
    out = str(tmp_path / "o.bcf")
    F.maf_filter(v, out, maf_min=0.20, stat_samples=",".join(cohort))
    assert _kept(out) == [1000]


def test_subset_and_per_group_frequencies_are_refused_together(tmp_path):
    v, cohort, _ = _vcf(tmp_path, [(1000, 3, 0)])
    meta = tmp_path / "m.tsv"
    meta.write_text("sample\tcountry\n" + "".join(f"{s}\tA\n" for s in cohort))
    with pytest.raises(SystemExit, match="one or the other"):
        F.maf_filter(v, str(tmp_path / "o.bcf"), meta=str(meta), group_col="country",
                     stat_samples=",".join(cohort))


# ---- locus_missingness_filter --------------------------------------------------------

def test_missingness_and_coverage_can_be_judged_on_the_cohort(tmp_path):
    """A locus the cohort covers well, inside a callset where the extras do not."""
    v, cohort, extra = _vcf(tmp_path, [(1000, 3, 0)])
    txt = pathlib.Path(v).read_text().splitlines()
    cols = txt[-1].split("\t")
    first = 9 + len(cohort) + 5                          # 9 fixed VCF columns, then samples
    for i in range(first, len(cols)):                    # most extras uncalled and at 0x
        cols[i] = "./.:0,0"
    pathlib.Path(v).write_text("\n".join(txt[:-1] + ["\t".join(cols)]) + "\n")
    with_ads = str(tmp_path / "ads.bcf")
    F.singleton_add_ads(v, with_ads, min_samples=0)
    plain = str(tmp_path / "plain.bcf")
    F.locus_missingness_filter(with_ads, plain)
    assert _kept(plain) == []                            # the extras sink it
    sub = str(tmp_path / "sub.bcf")
    F.locus_missingness_filter(with_ads, sub, stat_samples=_list_file(tmp_path, cohort))
    assert _kept(sub) == [1000]
    assert len(_samples(sub)) == len(cohort) + len(extra)


# ---- singleton_filter_add_ads --------------------------------------------------------

def test_a_carrier_count_can_be_judged_on_the_cohort(tmp_path):
    """1000 is carried by two cohort samples, 2000 by one cohort sample and many extras.
    'Carried by more than one sample' is a question about the analysis cohort."""
    v, cohort, _ = _vcf(tmp_path, [(1000, 2, 0), (2000, 1, 20)])
    plain = str(tmp_path / "plain.bcf")
    F.singleton_add_ads(v, plain, min_samples=1)
    assert _kept(plain) == [1000, 2000]
    sub = str(tmp_path / "sub.bcf")
    F.singleton_add_ads(v, sub, min_samples=1, stat_samples=_list_file(tmp_path, cohort))
    assert _kept(sub) == [1000]
    assert len(_samples(sub)) == 40


# ---- the plumbing --------------------------------------------------------------------

def test_a_name_not_in_the_callset_is_reported_not_fatal(tmp_path, capsys):
    v, cohort, _ = _vcf(tmp_path, [(1000, 3, 0)])
    out = str(tmp_path / "o.bcf")
    F.maf_filter(v, out, maf_min=0.20, stat_samples=",".join(cohort + ["ghost1", "ghost2"]))
    assert _kept(out) == [1000]
    cap = capsys.readouterr()
    assert "2 named sample(s) are not in the callset" in cap.out + cap.err


def test_no_named_sample_in_the_callset_is_an_error(tmp_path):
    v, _, _ = _vcf(tmp_path, [(1000, 3, 0)])
    with pytest.raises(SystemExit, match="none of the"):
        F.maf_filter(v, str(tmp_path / "o.bcf"), maf_min=0.2, stat_samples="nobody1,nobody2")


def test_an_unindexed_input_still_works(tmp_path):
    """The step is run by hand on another tool's output as often as inside a pipeline, and
    `bcftools annotate` wants an index on both sides."""
    v, cohort, _ = _vcf(tmp_path, [(1000, 3, 0)])
    plain_bcf = str(tmp_path / "noidx.bcf")
    subprocess.run(["bcftools", "view", v, "-Ob", "-o", plain_bcf], check=True,
                   stderr=subprocess.DEVNULL)
    assert not pathlib.Path(plain_bcf + ".csi").exists()
    out = str(tmp_path / "o.bcf")
    F.maf_filter(plain_bcf, out, maf_min=0.20, stat_samples=_list_file(tmp_path, cohort))
    assert _kept(out) == [1000]


def test_the_pipeline_threads_stat_samples_through_the_cohort_statistic_steps(tmp_path):
    v, cohort, _ = _vcf(tmp_path, [(1000, 3, 0), (2000, 0, 12)])
    cfg = {"stat_samples": _list_file(tmp_path, cohort),
           "steps": [{"name": "singleton_filter_add_ads", "params": {"min_samples": 0}},
                     {"name": "maf_filter", "params": {"maf_min": 0.20}}]}
    tally = P.run_pipeline(v, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    assert tally[-1]["variants"] == 1
    import json
    used = json.loads((tmp_path / "run" / "config_used.json").read_text())
    for step in used["steps"]:
        assert step["params"]["stat_samples"] == cfg["stat_samples"]


def test_the_sample_judging_steps_do_not_take_it():
    """`sample_coverage_filter` and `fws_filter` judge samples, not loci, so "whose
    statistic decides" does not arise -- give them fewer samples instead."""
    assert P.STAT_SUBSETTABLE == {"singleton_filter_add_ads", "locus_missingness_filter",
                                  "maf_filter"}
    for name in ("sample_coverage_filter", "fws_filter"):
        assert name not in P.STAT_SUBSETTABLE
    P.validate_config(P.DEFAULT_CONFIG)
    assert P.DEFAULT_CONFIG["stat_samples"] is None
