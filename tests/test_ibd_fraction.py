"""Per-pair IBD fraction: the edge-list columns and the fraction itself."""

from __future__ import annotations

import pandas as pd
import pytest

from plasgenomicsutils.lib import ibd_fraction as F


def _blocks(tmp_path, rows):
    p = tmp_path / "blocks.tsv"
    pd.DataFrame(rows, columns=["sample1", "sample2", "chr", "start", "end",
                                "different", "Nsnp"]).to_csv(p, sep="\t", index=False)
    return str(p)


def test_the_pair_table_is_an_edge_list(tmp_path):
    """sample1/sample2 are columns, not something to recover by splitting `pair` -- which is
    guesswork as soon as a sample name contains the separator."""
    b = _blocks(tmp_path, [
        ("a", "b", "c1", 0, 1_000_000, 0, 100),
        ("b", "a", "c1", 2_000_000, 2_500_000, 0, 100),   # same pair, listed the other way
        ("a", "c", "c1", 0, 100_000, 1, 100),             # compared, not IBD
    ])
    out = F.per_pair_fraction(b, sep="\t", bp_per_cm=15_000, callable_cm=100,
                              min_block_snp=0, min_block_kb=0)

    assert {"pair", "sample1", "sample2"} <= set(out.columns)
    assert (out["sample1"] < out["sample2"]).all()        # one orientation per pair
    assert (out["pair"] == out["sample1"] + "__" + out["sample2"]).all()
    assert len(out) == 2                                   # a-b and a-c, each once

    ab = out[out["pair"] == "a__b"].iloc[0]
    assert ab["sample1"] == "a" and ab["sample2"] == "b"
    # both segments are summed whichever way round they were listed. hmmibd-rs reports the
    # first and last SNP inclusive, and read_blocks shifts `end` to make the interval
    # half-open, so each segment is one base longer than end - start.
    assert ab["total_ibd_bp"] == (1_000_000 + 1) + (500_000 + 1)


def test_the_fraction_columns_say_what_their_denominator_is(tmp_path):
    b = _blocks(tmp_path, [("a", "b", "c1", 0, 1_500_000, 0, 100)])
    out = F.per_pair_fraction(b, sep="\t", bp_per_cm=15_000, callable_cm=100,
                              min_block_snp=0, min_block_kb=0)

    assert "ibd_fraction_accessible" in out.columns
    assert "f" not in out.columns                          # the old name is gone
    ab = out.iloc[0]
    assert ab["total_ibd_bp"] == 1_500_001                   # half-open: last SNP included
    assert ab["total_ibd_cm"] == pytest.approx(1_500_001 / 15_000)
    # the fraction is that map length over the callable one it is named for
    assert ab["ibd_fraction_accessible"] == pytest.approx(ab["total_ibd_cm"] / 100)


def test_a_compared_pair_with_no_IBD_still_gets_a_row(tmp_path):
    """The denominator of every downstream summary is the pairs that were compared, so a
    zero-sharing pair has to be present rather than absent."""
    b = _blocks(tmp_path, [
        ("a", "b", "c1", 0, 900_000, 0, 100),
        ("a", "c", "c1", 0, 900_000, 1, 100),      # different == 1: compared, not IBD
    ])
    out = F.per_pair_fraction(b, sep="\t", bp_per_cm=15_000, callable_cm=100,
                              min_block_snp=0, min_block_kb=0)
    ac = out[out["pair"] == "a__c"].iloc[0]
    assert ac["total_ibd_bp"] == 0
    assert ac["ibd_fraction_accessible"] == 0
    assert pd.isna(ac["gen_to_mrca_approx"])       # undefined with no shared segment


def test_short_segments_are_filtered_from_the_numerator_only(tmp_path):
    b = _blocks(tmp_path, [
        ("a", "b", "c1", 0, 5_000, 0, 2),          # too short and too few SNPs
    ])
    out = F.per_pair_fraction(b, sep="\t", bp_per_cm=15_000, callable_cm=100,
                              min_block_snp=15, min_block_kb=15)
    assert len(out) == 1                            # the pair survives
    assert out.iloc[0]["total_ibd_bp"] == 0         # its spurious segment does not
