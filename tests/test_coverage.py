"""Coverage statistics, checked against BAMs whose depth is known by construction."""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pysam
import pytest

from plasgenomicsutils.lib.coverage import (
    ChromDepth,
    annotate_regions,
    dropout_regions,
    load_bed,
    sample_coverage,
    summarise_histogram,
)

CHROM = "Pf3D7_01_v3"
CHROM_LEN = 2000
READ_LEN = 50


def _write_bam(path, intervals, sample="s1", chrom_len=CHROM_LEN, mapq=60):
    """One BAM whose depth is exactly `n` over each (start, end, n) interval.

    Reads are laid down one per base position so a stack of `n` reads covering
    [start, end) makes the depth there exactly n -- no need to reason about overlaps.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"LN": chrom_len, "SN": CHROM}],
              "RG": [{"ID": "rg1", "SM": sample}]}
    reads = []
    for start, end, n in intervals:
        for _ in range(n):
            for pos in range(start, end, READ_LEN):
                span = min(READ_LEN, end - pos)
                reads.append((pos, span))
    reads.sort()
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (pos, span) in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * span
            a.query_qualities = pysam.qualitystring_to_array("I" * span)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigarstring = f"{span}M"
            out.write(a)
    pysam.index(str(path))
    return str(path)


@pytest.fixture
def flat_bam(tmp_path):
    # depth 10 over [0, 1000), depth 0 over [1000, 2000)
    return _write_bam(tmp_path / "flat.bam", [(0, 1000, 10)])


def test_histogram_summary_matches_numpy():
    depth = np.array([0] * 10 + [4] * 30 + [9] * 60)
    hist = np.bincount(depth)
    s = summarise_histogram(hist, thresholds=(1, 5, 10))
    assert s["bases"] == 100
    assert s["mean"] == pytest.approx(depth.mean())
    assert s["sd"] == pytest.approx(depth.std())      # population SD
    assert s["median"] == 9
    assert s["min"] == 0 and s["max"] == 9
    assert s["pct_zero"] == pytest.approx(10.0)
    assert s["pct_ge_1x"] == pytest.approx(90.0)      # the 4s and 9s
    assert s["pct_ge_5x"] == pytest.approx(60.0)      # only the 9s
    assert s["pct_ge_10x"] == pytest.approx(0.0)


def test_empty_histogram_is_all_na():
    s = summarise_histogram(np.zeros(1, dtype=np.int64))
    assert s["bases"] == 0
    assert np.isnan(s["mean"]) and np.isnan(s["median"])


def test_depth_is_exact_over_a_known_bam(flat_bam):
    per_chrom, windows = sample_coverage(flat_bam, window=500, engine="pysam",
                                         thresholds=(1, 5, 10, 20))
    row = next(r for r in per_chrom if r["chrom"] == "genome")
    assert row["bases"] == CHROM_LEN
    assert row["mean"] == pytest.approx(10 * 1000 / CHROM_LEN)   # half the genome at 10x
    assert row["max"] == 10
    assert row["pct_zero"] == pytest.approx(50.0)
    assert row["pct_ge_10x"] == pytest.approx(50.0)
    assert row["pct_ge_20x"] == pytest.approx(0.0)
    # the covered half splits into two full windows at 10x, the empty half into two at 0
    w = pd.DataFrame(windows).sort_values("start")
    assert list(w["mean_depth"]) == [10.0, 10.0, 0.0, 0.0]
    assert set(w["bases"]) == {500}


def test_regions_restrict_the_denominator(tmp_path, flat_bam):
    bed = tmp_path / "core.bed"
    bed.write_text(f"{CHROM}\t0\t1000\n")          # only the covered half
    per_chrom, windows = sample_coverage(flat_bam, regions=load_bed(bed), window=500,
                                         engine="pysam")
    row = next(r for r in per_chrom if r["chrom"] == "genome")
    assert row["bases"] == 1000                     # bases outside the BED are not counted
    assert row["mean"] == pytest.approx(10.0)       # ...so the mean is the in-region mean
    assert row["pct_zero"] == pytest.approx(0.0)
    assert len(windows) == 2                        # empty windows are dropped


def test_overlapping_bed_intervals_count_each_base_once(tmp_path, flat_bam):
    bed = tmp_path / "overlap.bed"
    bed.write_text(f"{CHROM}\t0\t600\n{CHROM}\t400\t1000\n")
    merged = load_bed(bed)
    assert merged[CHROM] == [(0, 1000)]
    row = next(r for r in sample_coverage(flat_bam, regions=merged, engine="pysam")[0]
               if r["chrom"] == "genome")
    assert row["bases"] == 1000


def test_min_mapq_drops_low_quality_reads(tmp_path):
    bam = _write_bam(tmp_path / "lowq.bam", [(0, 1000, 10)], mapq=5)
    keep = next(r for r in sample_coverage(bam, engine="pysam", min_mapq=1)[0]
                if r["chrom"] == "genome")
    drop = next(r for r in sample_coverage(bam, engine="pysam", min_mapq=30)[0]
                if r["chrom"] == "genome")
    assert keep["max"] == 10
    assert drop["max"] == 0          # every read is below the floor


def test_chunking_does_not_change_the_answer(flat_bam):
    whole = sample_coverage(flat_bam, engine="pysam", chunk=10_000)[0]
    split = sample_coverage(flat_bam, engine="pysam", chunk=317)[0]   # ragged on purpose
    assert whole == split


def test_chrom_depth_accumulates_across_chunks():
    cd = ChromDepth(CHROM, 100, window=50)
    cd.add(np.full(50, 3, dtype=np.int64), 0)
    cd.add(np.full(50, 7, dtype=np.int64), 50)
    assert cd.hist[3] == 50 and cd.hist[7] == 50
    assert list(cd.win_depth_sum / cd.win_bases) == [3.0, 7.0]


# --------------------------------------------------------------------------- #
#  Dropouts                                                                    #
# --------------------------------------------------------------------------- #


def _windows(spec):
    """spec: {sample: [(start, mean_depth), ...]} on a 100 bp grid."""
    rows = []
    for sample, wins in spec.items():
        for start, depth in wins:
            rows.append({"sample": sample, "chrom": CHROM, "start": start,
                         "end": start + 100, "bases": 100, "mean_depth": depth})
    return pd.DataFrame(rows)


def test_dropout_needs_almost_every_sample_to_be_uncovered():
    # window 0: everyone dead. window 100: one sample fine, so not a dropout.
    w = _windows({
        "a": [(0, 0.0), (100, 0.0)],
        "b": [(0, 0.1), (100, 0.0)],
        "c": [(0, 0.0), (100, 40.0)],
    })
    out = dropout_regions(w, min_depth=5, min_frac_samples=0.9)
    assert list(out["start"]) == [0]
    assert out.loc[0, "frac_samples_uncovered"] == pytest.approx(1.0)
    assert out.loc[0, "n_samples"] == 3

    # relax the fraction and the second window qualifies too -- and being adjacent, the
    # two come back as one region rather than two
    out2 = dropout_regions(w, min_depth=5, min_frac_samples=0.6)
    assert list(zip(out2["start"], out2["end"])) == [(0, 200)]
    assert out2.loc[0, "n_windows"] == 2


def test_adjacent_dropout_windows_merge_into_one_region():
    w = _windows({s: [(0, 0.0), (100, 0.0), (200, 50.0), (300, 0.0)]
                  for s in ("a", "b")})
    out = dropout_regions(w, min_depth=5, min_frac_samples=0.9).sort_values("start")
    assert list(zip(out["start"], out["end"])) == [(0, 200), (300, 400)]
    assert list(out["n_windows"]) == [2, 1]
    assert list(out["length"]) == [200, 100]

    # a gap of one covered window can be bridged, and short regions dropped
    bridged = dropout_regions(w, min_depth=5, min_frac_samples=0.9, merge_gap=100)
    assert list(zip(bridged["start"], bridged["end"])) == [(0, 400)]
    assert dropout_regions(w, min_depth=5, min_frac_samples=0.9,
                           min_length=150)["length"].tolist() == [200]


def test_no_dropouts_returns_an_empty_frame():
    w = _windows({"a": [(0, 30.0)], "b": [(0, 25.0)]})
    out = dropout_regions(w, min_depth=5, min_frac_samples=0.9)
    assert out.empty


def test_regions_are_annotated_with_overlapping_genes():
    regions = pd.DataFrame([{"chrom": CHROM, "start": 100, "end": 400}])
    genes = pd.DataFrame([
        {"name": "geneA", "chr": CHROM, "start": 0, "end": 150},     # overlaps the start
        {"name": "geneB", "chr": CHROM, "start": 380, "end": 900},   # overlaps the end
        {"name": "geneC", "chr": CHROM, "start": 400, "end": 500},   # abuts, half-open
    ])
    out = annotate_regions(regions, genes)
    assert out.loc[0, "genes"] == "geneA,geneB"


# --------------------------------------------------------------------------- #
#  Engine agreement                                                            #
# --------------------------------------------------------------------------- #


def test_the_engines_agree_when_no_mates_overlap(flat_bam):
    """The engines define depth differently -- mosdepth counts a fragment once where the
    mates overlap, pysam counts both reads (matching `samtools depth`). On real WGS that
    puts mosdepth 2-3% low. These fixture reads are unpaired, so nothing overlaps and the
    two must agree exactly; anywhere else they must not be assumed interchangeable."""
    from plasgenomicsutils.lib.coverage import mosdepth_available

    if not mosdepth_available():
        pytest.skip("mosdepth not installed")
    a = sample_coverage(flat_bam, window=500, engine="pysam")
    b = sample_coverage(flat_bam, window=500, engine="mosdepth")
    strip = lambda rows: [{k: v for k, v in r.items() if k != "engine"} for r in rows]
    assert strip(a[0]) == strip(b[0])
    assert a[1] == b[1]


def test_the_engine_used_is_recorded_in_every_row(flat_bam):
    """Two engines that count different things must never be silently interchangeable."""
    rows, _ = sample_coverage(flat_bam, engine="pysam")
    assert {r["engine"] for r in rows} == {"pysam"}


def test_per_base_rle_expands_to_the_same_depths(tmp_path):
    """The mosdepth reader is exercised even without the binary, via a hand-written file."""
    from plasgenomicsutils.lib.coverage import _depth_from_per_base

    p = tmp_path / "s.per-base.bed.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(f"{CHROM}\t0\t10\t3\n{CHROM}\t10\t25\t0\n{CHROM}\t25\t30\t7\n")
    depth = _depth_from_per_base(p, {CHROM: 30})[CHROM]
    assert depth.tolist() == [3] * 10 + [0] * 15 + [7] * 5
