"""The two end-of-chain reports: what the callset looks like, with nothing changed."""

import pathlib
import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import callset_summary as CS
from plasgenomicsutils.lib import filter_pipeline as P
from plasgenomicsutils.lib import vcf_filters as F

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not on PATH")

DATA = pathlib.Path(__file__).parent / "data"
BCF = DATA / "ghana_cambodia.pf7.tiny.bcf"


def _vcf(tmp_path, records):
    """(pos, ref, alt) records, one sample, so the class and allele count are the whole story."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"]
    for pos, ref, alt in records:
        hdr.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\t.\t.\tGT\t0/1")
    p = tmp_path / "v.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def test_variant_summary_counts_by_class_and_allele_number_with_both_fractions(tmp_path):
    v = _vcf(tmp_path, [(1, "A", "T"), (2, "A", "T"), (3, "A", "T,G"), (4, "AT", "A"),
                        (5, "A", "AT,ATT,ATTT")])
    rows = {(r["class"], r["n_alt"]): r for r in CS.variant_summary_table(v)}
    assert rows[("snps", "all")]["count"] == 3 and rows[("snps", "all")]["frac_total"] == 0.6
    assert rows[("snps", 1)] == {"class": "snps", "n_alt": 1, "alleles": "biallelic",
                                 "count": 2, "frac_total": 0.4, "frac_class": 0.6667}
    assert rows[("snps", 2)]["alleles"] == "triallelic" and rows[("snps", 2)]["count"] == 1
    assert rows[("indels", 3)]["alleles"] == "tetra-allelic"
    assert rows[("total", "all")]["count"] == 5
    note = CS.variant_summary_note(CS.variant_summary_table(v))
    assert note.startswith("5 record(s): snps 3 60.0% (biallelic 2 66.7%, triallelic 1 33.3%)")


def test_variant_summary_on_an_empty_callset_is_one_zero_row(tmp_path):
    v = _vcf(tmp_path, [])
    rows = CS.variant_summary_table(v)
    assert rows == [{"class": "total", "n_alt": "all", "alleles": "all", "count": 0,
                     "frac_total": 0.0, "frac_class": 0.0}]


def test_sample_summary_pairs_coverage_with_fws_and_drops_nobody(tmp_path):
    rows, n_sites = CS.sample_summary_table(str(BCF))      # no ADS: added on a temp copy
    assert len(rows) == 60 and n_sites > 0
    assert {r["sample"] for r in rows} == set(
        subprocess.run(["bcftools", "query", "-l", str(BCF)], capture_output=True,
                       text=True).stdout.split())
    r = rows[0]
    assert set(r) == set(CS.SAMPLE_SUMMARY_COLUMNS)
    assert 0 <= r["frac_covered"] <= 1
    assert r["fws"] is None or 0 <= r["fws"] <= 1
    # the columns agree with the filters they mirror, at the same thresholds
    cov = {c["sample"]: c["dropped"] for c in F.sample_coverage_table(
        str(tmp_path / "ads.bcf") if F.singleton_add_ads(str(BCF), str(tmp_path / "ads.bcf")) is not None else "")}
    assert all(r["would_drop_coverage"] == cov[r["sample"]] for r in rows)
    out = tmp_path / "s.tsv"
    CS.write_sample_summary(rows, str(out))
    assert out.read_text().splitlines()[0] == "\t".join(CS.SAMPLE_SUMMARY_COLUMNS)


def test_the_default_chain_ends_with_both_reports_and_they_run(tmp_path):
    names = [s["name"] for s in P.DEFAULT_CONFIG["steps"]]
    assert names[-2:] == ["sample_summary", "variant_summary"]
    assert all(s.get("report") for s in P.DEFAULT_CONFIG["steps"][-2:])
    P.validate_config(P.DEFAULT_CONFIG)
    cfg = {"steps": [{"name": "singleton_filter_add_ads"},
                     {"name": "sample_summary", "report": True, "ext": "tsv",
                      "params": {"frac_min": 0.5}},
                     {"name": "variant_summary", "report": True, "ext": "tsv"}]}
    tally = P.run_pipeline(str(BCF), str(tmp_path), cfg, emit_snp_bed=False)
    reports = [t for t in tally if t.get("report")]
    assert [r["step"] for r in reports] == ["sample_summary", "variant_summary"]
    assert reports[0]["rows"] == 60
    assert pathlib.Path(reports[1]["path"]).read_text().startswith("class\tn_alt")
    # a report leaves the callset alone: the filter output is the last callset
    assert tally[-1]["step"] == "variant_summary" and "variants" not in tally[-1]


def test_a_report_param_is_validated_against_the_table_builder(tmp_path):
    with pytest.raises(SystemExit, match="sample_summary"):
        P.validate_config({"steps": [{"name": "sample_summary", "report": True,
                                      "params": {"frac_mn": 0.5}}]})


# ---- where the variation sits on the frequency spectrum -------------------------------

def _freq_vcf(tmp_path, specs, n=100):
    """`specs` are (pos, n_carriers_groupA, n_carriers_groupB) with 50 samples per group."""
    half = n // 2
    a = [f"a{i:02d}" for i in range(half)]
    b = [f"b{i:02d}" for i in range(half)]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(a + b)]
    for pos, na, nb in specs:
        cols = ["1/1" if i < na else "0/0" for i in range(half)]
        cols += ["1/1" if i < nb else "0/0" for i in range(half)]
        hdr.append(f"chr1\t{pos}\t.\tA\tT\t500\t.\t.\tGT\t" + "\t".join(cols))
    p = tmp_path / "freq.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p), a, b


def test_the_spectrum_counts_records_at_each_mark(tmp_path):
    # 50 samples per group, 100 total: one carrier is 1% of the 200 alleles
    v, _, _ = _freq_vcf(tmp_path, [(1000, 1, 0), (2000, 2, 0), (3000, 5, 0), (4000, 20, 0)])
    rows = CS.maf_spectrum_table(v)
    at = {r["maf_min"]: r["at_or_above"] for r in rows}
    assert at[0.01] == 4 and at[0.02] == 3 and at[0.05] == 2 and at[0.10] == 1
    band = {r["maf_min"]: r["in_band_below_next"] for r in rows}
    assert band[0.01] == 1 and band[0.02] == 1 and band[0.05] == 1 and band[0.10] == 1
    assert all(r["group"] == "ALL" and r["n_samples"] == 100 for r in rows)


def test_the_spectrum_is_computed_per_group_when_asked(tmp_path):
    """A grouped floor keeps a record when any one group clears it, so the per-group
    spectrum is the one that predicts what a grouped filter does."""
    # hom-alt carriers contribute two alleles: 5 carriers of 50 group-A samples is 10 of
    # that group's 100 alleles (10%), and 10 of the cohort's 200 (5%)
    v, a, b = _freq_vcf(tmp_path, [(1000, 5, 0)])
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tcountry\n" + "".join(f"{s}\tA\n" for s in a)
                    + "".join(f"{s}\tB\n" for s in b))
    rows = CS.maf_spectrum_table(v, meta=str(meta), group_col="country")
    by = {(r["group"], r["maf_min"]): r["at_or_above"] for r in rows}
    assert by[("ALL", 0.05)] == 1 and by[("ALL", 0.10)] == 0     # 5% pooled
    assert by[("A", 0.10)] == 1                                  # 10% within its own group
    assert by[("B", 0.01)] == 0                                  # absent from the other
    assert {r["group"] for r in rows} == {"ALL", "A", "B"}
    note = CS.maf_spectrum_note(rows)
    assert note.startswith("ALL (n=100") and "A (n=50" in note


def test_the_spectrum_runs_before_the_frequency_filter_in_the_default_chain():
    names = [s["name"] for s in P.DEFAULT_CONFIG["steps"]]
    assert names.index("maf_spectrum") < names.index("maf_filter")
    assert P.DEFAULT_CONFIG["steps"][names.index("maf_spectrum")].get("report") is True


def test_the_default_frequency_floor_is_one_percent():
    """The stored callset uses the most permissive floor any downstream analysis wants: a
    stricter filter is one command away on a 1% callset and impossible on a 2% one."""
    step = next(s for s in P.DEFAULT_CONFIG["steps"] if s["name"] == "maf_filter")
    assert step["params"]["maf_min"] == 0.01
