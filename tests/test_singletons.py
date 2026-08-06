"""Per-sample singleton counts and the MAD outlier flag."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from plasgenomicsutils.lib.singletons import count_singletons, flag_outliers

DATA = Path(__file__).parent / "data"
BCF = DATA / "ghana_cambodia.pf7.tiny.bcf"


def _write_vcf(path, samples, rows):
    """rows: list of genotype strings per sample, e.g. ['0', '1', '.']."""
    head = ["##fileformat=VCFv4.2", '##contig=<ID=chr1,length=100000>',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    lines = list(head)
    for i, gts in enumerate(rows, start=1):
        lines.append(f"chr1\t{i * 100}\t.\tA\tT\t100\tPASS\t.\tGT\t" + "\t".join(gts))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_only_a_lone_carrier_counts_as_a_singleton(tmp_path):
    samples = ["a", "b", "c", "d"]
    vcf = _write_vcf(tmp_path / "t.vcf", samples, [
        ["1", "0", "0", "0"],      # singleton for a
        ["0", "1", "0", "0"],      # singleton for b
        ["1", "1", "0", "0"],      # two carriers: not a singleton
        ["1", "1", "1", "1"],      # fixed: not a singleton
        ["0", "0", "0", "0"],      # invariant
        ["1", "0", "0", "0"],      # another for a
    ])
    df, n_variants = count_singletons(vcf)
    assert n_variants == 6
    got = dict(zip(df["sample"], df["n_singleton"]))
    assert got == {"a": 2, "b": 1, "c": 0, "d": 0}
    assert list(df["n_called"]) == [6] * 4
    assert df.loc[df["sample"] == "a", "singleton_rate"].iloc[0] == pytest.approx(2000 / 6)


def test_missing_genotypes_do_not_count_as_reference(tmp_path):
    samples = ["a", "b", "c"]
    vcf = _write_vcf(tmp_path / "t.vcf", samples, [
        ["1", ".", "."],           # a is the only *carrier*, but also the only call
        ["1", "0", "."],
    ])
    df, _ = count_singletons(vcf, max_missing_frac=1.0)
    called = dict(zip(df["sample"], df["n_called"]))
    assert called == {"a": 2, "b": 1, "c": 0}      # missing is missing, not hom-ref


def test_thin_variants_are_skipped(tmp_path):
    samples = ["a", "b", "c", "d"]
    vcf = _write_vcf(tmp_path / "t.vcf", samples, [
        ["1", ".", ".", "."],      # 75% missing -- a private allele here means nothing
        ["1", "0", "0", "0"],
    ])
    lenient, n_all = count_singletons(vcf, max_missing_frac=1.0)
    strict, n_strict = count_singletons(vcf, max_missing_frac=0.5)
    assert n_all == 2 and n_strict == 1
    assert lenient.loc[lenient["sample"] == "a", "n_singleton"].iloc[0] == 2
    assert strict.loc[strict["sample"] == "a", "n_singleton"].iloc[0] == 1


def test_singleton_status_is_relative_to_the_samples_analysed(tmp_path):
    samples = ["a", "b", "c"]
    vcf = _write_vcf(tmp_path / "t.vcf", samples, [["1", "1", "0"]])
    full, _ = count_singletons(vcf)
    assert full["n_singleton"].sum() == 0          # two carriers in the full cohort
    subset, _ = count_singletons(vcf, samples=["a", "c"])
    assert subset.loc[subset["sample"] == "a", "n_singleton"].iloc[0] == 1


def test_outlier_flag_uses_the_median_not_the_mean():
    # one wild sample; a mean/SD rule would be dragged towards it, a MAD rule is not
    df = pd.DataFrame({"sample": list("abcdefghij"),
                       "singleton_rate": [1.0, 1.1, 0.9, 1.0, 1.2, 0.8, 1.0, 1.1, 0.9, 50.0]})
    out = flag_outliers(df, mad_cutoff=5)
    assert out.loc[out["sample"] == "j", "outlier"].iloc[0]
    assert not out.loc[out["sample"] != "j", "outlier"].any()


def test_a_uniform_cohort_flags_nobody():
    df = pd.DataFrame({"sample": list("abcd"), "singleton_rate": [2.0] * 4})
    out = flag_outliers(df)
    assert not out["outlier"].any()
    assert (out["mad_score"] == 0).all()           # zero MAD must not divide by zero


def test_runs_on_the_public_fixture():
    df, n_variants = count_singletons(str(BCF), max_missing_frac=0.2)
    assert n_variants > 0
    assert len(df) == 60
    assert (df["n_singleton"] >= 0).all()
    assert (df["n_singleton"] <= n_variants).all()
    flagged = flag_outliers(df)
    assert set(flagged.columns) >= {"mad_score", "outlier"}


def _write_vcf_ad(path, samples, rows):
    """rows: list of (gt, ad_depth) per sample, so depth-based calling can be tested."""
    head = ["##fileformat=VCFv4.2", '##contig=<ID=chr1,length=100000>',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    lines = list(head)
    for i, cells in enumerate(rows, start=1):
        fmt = [f"{gt}:{d},0" if gt == "0" else f"{gt}:0,{d}" for gt, d in cells]
        lines.append(f"chr1\t{i * 100}\t.\tA\tT\t100\tPASS\t.\tGT:AD\t" + "\t".join(fmt))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_a_zero_depth_reference_call_is_not_a_confident_call(tmp_path):
    """gVCF-derived callsets emit 0/0 with no reads rather than ./. -- a real pattern
    that would otherwise be counted as a called reference genotype."""
    samples = ["a", "b", "c"]
    vcf = _write_vcf_ad(tmp_path / "t.vcf", samples, [
        [("1", 30), ("0", 0), ("0", 25)],      # b has no reads at all
        [("0", 30), ("0", 0), ("0", 25)],
    ])
    lenient, _ = count_singletons(vcf, min_depth=0)
    strict, _ = count_singletons(vcf, min_depth=5)
    assert dict(zip(lenient["sample"], lenient["n_called"])) == {"a": 2, "b": 2, "c": 2}
    assert dict(zip(strict["sample"], strict["n_called"])) == {"a": 2, "b": 0, "c": 2}
    assert strict.attrs["n_low_depth"] == 2


def test_a_near_identical_pair_is_found_from_its_doubletons(tmp_path):
    # x and y carry the same 20 rare variants; nobody else carries any of them
    samples = ["x", "y"] + [f"o{i}" for i in range(6)]
    rows = []
    for _ in range(20):
        rows.append([("1", 30), ("1", 30)] + [("0", 30)] * 6)
    for i in range(6):                                   # each other sample has its own
        cells = [("0", 30)] * 8
        cells[2 + i] = ("1", 30)
        rows.append(cells)
    vcf = _write_vcf_ad(tmp_path / "dup.vcf", samples, rows)
    df, _ = count_singletons(vcf, min_depth=5)
    df = flag_outliers(df)

    row = df[df["sample"] == "x"].iloc[0]
    assert row["n_singleton"] == 0                       # its variants are never private
    assert row["top_partner"] == "y"
    assert row["frac_doubletons_with_partner"] == 1.0
    assert "near-identical to y" in row["flag"]
    assert df[df["sample"] == "y"].iloc[0]["flag"].endswith("near-identical to x")
    # the unrelated samples keep their singletons and are not flagged
    others = df[df["sample"].str.startswith("o")]
    assert (others["n_singleton"] == 1).all()
    assert (others["flag"] == "").all()


def test_merely_related_samples_are_not_called_near_identical():
    # a pair sharing 70% of its doubletons is related, not a duplicate -- the real
    # cohort has pairs at 0.58-0.78 that must not be flagged
    df = pd.DataFrame({
        "sample": ["a", "b", "c"], "singleton_rate": [1.0, 1.0, 1.0],
        "top_partner": ["b", "a", "a"],
        "frac_doubletons_with_partner": [0.70, 0.78, 0.05]})
    out = flag_outliers(df)
    assert (out["flag"] == "").all()
    # ...but a 0.95 pair is
    df.loc[0, "frac_doubletons_with_partner"] = 0.95
    assert "near-identical to b" in flag_outliers(df).iloc[0]["flag"]
