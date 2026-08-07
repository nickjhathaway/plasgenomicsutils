"""Pair-level IBD adjacency list over genes."""

import numpy as np
import pandas as pd

from plasgenomicsutils.lib.ibd_gene_pairs import (
    COLUMNS,
    gene_ibd_pairs as _gene_ibd_pairs,
)
# These fixtures use deliberately tiny blocks to make the overlap arithmetic checkable by
# hand, so they opt out of the default short-segment filter; the filter has its own tests.


def gene_ibd_pairs(*args, **kw):
    kw.setdefault("min_block_snp", 0)
    kw.setdefault("min_block_kb", 0)
    return _gene_ibd_pairs(*args, **kw)



GENES = pd.DataFrame({
    "name": ["g1", "g2"],
    "chr": ["Pf3D7_07_v3", "Pf3D7_07_v3"],
    "start": [1000, 5000],
    "end": [2000, 6000],
    "gene_id": ["PF3D7_g1", "PF3D7_g2"],
})


def _blocks(rows):
    return pd.DataFrame(rows, columns=["sample1", "sample2", "chr", "start", "end", "different"])


def test_complete_vs_partial_and_covered_span():
    b = _blocks([
        ("a", "b", "Pf3D7_07_v3", 500, 2500, 0),    # spans g1 entirely
        ("c", "d", "Pf3D7_07_v3", 1500, 2500, 0),   # covers the back half of g1
        ("e", "f", "Pf3D7_07_v3", 500, 1250, 0),    # covers the front quarter of g1
    ])
    out = gene_ibd_pairs(b, GENES.iloc[[0]])
    assert list(out.columns) == COLUMNS
    got = out.set_index("sample1")

    assert got.loc["a", "coverage"] == "complete"
    # a complete overlap reports the gene's own bounds, not the block's
    assert (got.loc["a", "covered_start"], got.loc["a", "covered_end"]) == (1000, 2000)
    assert got.loc["a", "percent_covered"] == 100.0

    assert got.loc["c", "coverage"] == "partial"
    assert (got.loc["c", "covered_start"], got.loc["c", "covered_end"]) == (1500, 2000)
    assert got.loc["c", "percent_covered"] == 50.0

    assert got.loc["e", "covered_bp"] == 250
    assert got.loc["e", "percent_covered"] == 25.0


def test_only_pairs_with_ibd_are_listed():
    b = _blocks([
        ("a", "b", "Pf3D7_07_v3", 500, 2500, 0),
        ("c", "d", "Pf3D7_07_v3", 500, 2500, 1),      # different==1 -> not IBD
        ("e", "f", "Pf3D7_07_v3", 3000, 4000, 0),     # IBD but nowhere near a gene
    ])
    out = gene_ibd_pairs(b, GENES)
    assert out["sample1"].tolist() == ["a"]


def test_several_genes_come_back_in_one_table():
    b = _blocks([
        ("a", "b", "Pf3D7_07_v3", 500, 6500, 0),      # spans both genes
        ("c", "d", "Pf3D7_07_v3", 5500, 6500, 0),     # only g2, partially
    ])
    out = gene_ibd_pairs(b, GENES)
    assert sorted(out["gene"].unique()) == ["g1", "g2"]
    assert len(out) == 3
    assert out[out["gene"] == "g2"]["sample1"].tolist() == ["a", "c"]
    assert out["gene_id"].tolist() == ["PF3D7_g1", "PF3D7_g2", "PF3D7_g2"]


def test_pairs_are_order_normalised_and_a_pair_may_repeat():
    b = _blocks([
        ("z", "a", "Pf3D7_07_v3", 500, 1200, 0),      # reported z,a -> stored a,z
        ("a", "z", "Pf3D7_07_v3", 1800, 2500, 0),     # a second segment over the same gene
    ])
    out = gene_ibd_pairs(b, GENES.iloc[[0]])
    assert out["sample1"].tolist() == ["a", "a"]
    assert out["sample2"].tolist() == ["z", "z"]
    assert sorted(out["covered_bp"]) == [200, 200]


def test_within_pads_selection_but_not_coverage():
    b = _blocks([("a", "b", "Pf3D7_07_v3", 2100, 2200, 0)])   # just past g1's end
    assert gene_ibd_pairs(b, GENES.iloc[[0]]).empty
    out = gene_ibd_pairs(b, GENES.iloc[[0]], within=500)
    assert len(out) == 1
    assert out["covered_bp"].iloc[0] == 0                     # nothing of the gene itself
    assert out["percent_covered"].iloc[0] == 0.0
    assert np.isnan(out["covered_start"].iloc[0])


def test_no_overlap_returns_an_empty_frame_with_the_columns():
    b = _blocks([("a", "b", "Pf3D7_08_v3", 500, 2500, 0)])    # different chromosome
    out = gene_ibd_pairs(b, GENES)
    assert out.empty and list(out.columns) == COLUMNS
