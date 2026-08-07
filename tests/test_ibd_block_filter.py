"""Short / SNP-poor IBD segments are dropped by default (commonly spurious)."""

import pandas as pd
import pytest

from plasgenomicsutils.lib.ibd_matrix import (
    IBD_MIN_BLOCK_KB,
    IBD_MIN_BLOCK_SNP,
    filter_ibd_blocks,
)


def _blocks(rows):
    return pd.DataFrame(rows, columns=["sample1", "sample2", "chr", "start", "end", "Nsnp"])


def test_defaults_are_the_conventional_thresholds():
    assert (IBD_MIN_BLOCK_SNP, IBD_MIN_BLOCK_KB) == (15, 15.0)


def test_drops_short_or_snp_poor_segments():
    b = _blocks([
        ("a", "b", "c1", 0, 20_000, 30),      # long enough, SNP-rich  -> keep
        ("c", "d", "c1", 0, 14_999, 30),      # 1 bp too short         -> drop
        ("e", "f", "c1", 0, 20_000, 14),      # one SNP too few        -> drop
        ("g", "h", "c1", 0, 15_000, 15),      # exactly at both bounds -> keep
    ])
    kept = filter_ibd_blocks(b)
    assert kept["sample1"].tolist() == ["a", "g"]


def test_either_criterion_can_be_disabled():
    b = _blocks([("a", "b", "c1", 0, 1_000, 30), ("c", "d", "c1", 0, 20_000, 2)])
    assert filter_ibd_blocks(b, min_kb=0)["sample1"].tolist() == ["a"]      # SNP filter only
    assert filter_ibd_blocks(b, min_snp=0)["sample1"].tolist() == ["c"]     # length filter only
    assert len(filter_ibd_blocks(b, min_snp=0, min_kb=0)) == 2              # both off


def test_missing_Nsnp_warns_and_still_filters_on_length():
    b = _blocks([("a", "b", "c1", 0, 20_000, 30)]).drop(columns=["Nsnp"])
    with pytest.warns(UserWarning, match="no 'Nsnp' column"):
        kept = filter_ibd_blocks(b)
    assert len(kept) == 1


def test_filtering_does_not_shrink_the_compared_pair_set():
    """The denominator comes from the unfiltered frame, so a pair with only a short
    segment still counts as compared -- it just contributes no IBD."""
    from plasgenomicsutils.lib.ibd_gene_overlap import gene_block_overlap

    blocks = pd.DataFrame({
        "sample1": ["a", "c"], "sample2": ["b", "d"], "chr": ["Pf3D7_07_v3"] * 2,
        "start": [0, 0], "end": [40_000, 900], "different": [0, 0], "Nsnp": [40, 3],
    })
    genes = pd.DataFrame({"name": ["g"], "chr": ["Pf3D7_07_v3"],
                          "start": [100], "end": [800]})
    s2g = {"a": "A", "b": "A", "c": "B", "d": "B"}
    out = gene_block_overlap(blocks, genes, s2g)
    bb = out[(out["group_a"] == "B") & (out["group_b"] == "B")].iloc[0]
    assert bb["n_pairs_total"] == 1        # c,d still counted as a compared pair
    assert bb["n_pairs_ibd"] == 0          # but their short segment is not evidence
