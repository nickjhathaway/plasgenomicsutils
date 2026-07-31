"""Tests for the strand-bias artifact core (the BAM-free, pure-numeric paths)."""

import pytest

pytest.importorskip("scipy")

from plasgenomicsutils.lib.strandbias import (
    scan_reference_triggers, strand_bias_verdict, summarize_reverse_alt,
)


def test_verdict_flags_a_strand_restricted_artifact():
    # archetype: ~45% alt on reverse, ~0% on forward, deep on both strands -> drop
    v = strand_bias_verdict(ref_fwd=1500, alt_fwd=2, ref_rev=1200, alt_rev=1000)
    assert v["drop"] is True
    assert v["ratio"] < 0.15
    assert v["sb_phred"] > 60
    assert v["reasons"]


def test_verdict_keeps_a_strand_balanced_het():
    # equal alt fraction on both strands -> real het, no bias
    v = strand_bias_verdict(ref_fwd=50, alt_fwd=50, ref_rev=50, alt_rev=50)
    assert v["drop"] is False
    assert v["sb_phred"] == pytest.approx(0.0, abs=1e-9)
    assert v["ratio"] == pytest.approx(1.0)


def test_verdict_min_minor_depth_guard():
    # ratio is 0 (all alt on reverse) but the forward strand is too shallow to have
    # detected the alt; with the Fisher rule disabled, the depth guard must hold.
    kept = strand_bias_verdict(5, 0, 5, 45, sb_hard=1e9)
    assert kept["drop"] is False           # minor_depth 5 < default 20 -> not flagged
    dropped = strand_bias_verdict(5, 0, 5, 45, sb_hard=1e9, min_minor_depth=4)
    assert dropped["drop"] is True         # now 5 >= 4 -> strand-restricted fires


def test_verdict_low_bq_soft_rule():
    # softer strand skew alone does not drop; low alt BQ corroboration flips it
    kept = strand_bias_verdict(40, 2, 20, 8, sb_hard=1e9)     # ratio ~0.26, no BQ
    assert kept["drop"] is False
    dropped = strand_bias_verdict(40, 2, 20, 8, sb_hard=1e9, alt_bq_median=9)
    assert dropped["drop"] is True


def test_scan_reference_triggers_finds_motifs():
    # GGGG homopolymer (>=3), a GGC motif, and a GC-rich island in an AT background
    seq = "AAAAGGGGAAAAGGCAAAAGCGCGCGCGCGCGCGCGCGCAAAA"
    feats = scan_reference_triggers(seq, start=1000)
    kinds = {f["kind"] for f in feats}
    assert "homopolymer" in kinds
    assert "motif" in kinds
    assert "gc_island" in kinds
    homo = next(f for f in feats if f["kind"] == "homopolymer")
    assert homo["start"] == 1004 and homo["detail"] == "GGGG"   # 0-based genomic coord


def _rev_alt(cycle, bq=9, start=975431):
    return dict(strand="-", is_alt=True, base="A", bq=bq, seq_cycle=cycle,
                dist_nearest_end=30, read_len=150, sc_5p=0, sc_3p=0, has_sa=0,
                seq_start=start, seq_start_ext=start)


def test_summarize_reverse_alt_computes_strand_vaf():
    rows = (
        [_rev_alt(c) for c in (10, 40, 70, 100, 130)]                 # 5 reverse ALT
        + [dict(strand="-", is_alt=False, base="G", bq=35, seq_cycle=50,
                dist_nearest_end=40, read_len=150, sc_5p=0, sc_3p=0, has_sa=0,
                seq_start=975431, seq_start_ext=975431) for _ in range(5)]  # 5 reverse REF
        + [dict(strand="+", is_alt=False, base="G", bq=35, seq_cycle=50,
                dist_nearest_end=40, read_len=150, sc_5p=0, sc_3p=0, has_sa=0,
                seq_start=975200, seq_start_ext=975200) for _ in range(10)]  # 10 fwd REF
    )
    s = summarize_reverse_alt(rows, alt_base="A")
    assert s["n_rev_alt"] == 5 and s["n_fwd_alt"] == 0
    assert s["vaf_rev"] == pytest.approx(50.0)   # 5 alt / 10 reverse
    assert s["vaf_fwd"] == pytest.approx(0.0)
    assert s["rev_alt"]["bq_median"] == 9        # low BQ = the SSE signature
    assert s["rev_alt"]["modal_start_frac"] == pytest.approx(1.0)  # all share priming start


def test_summarize_handles_no_reverse_alt():
    rows = [dict(strand="+", is_alt=True, base="A", bq=35, seq_cycle=50,
                 dist_nearest_end=40, read_len=150, sc_5p=0, sc_3p=0, has_sa=0,
                 seq_start=1, seq_start_ext=1)]
    s = summarize_reverse_alt(rows)
    assert s["rev_alt"] is None
    assert s["n_fwd_alt"] == 1
