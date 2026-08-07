"""Per-gene IBD-block overlap: numerator = distinct pairs overlapping, denominator = all pairs."""

import numpy as np
import pandas as pd

from plasgenomicsutils.lib.ibd_gene_overlap import (
    gene_block_overlap as _gene_block_overlap,
)
# These fixtures use deliberately tiny blocks to make the overlap arithmetic checkable by
# hand, so they opt out of the default short-segment filter; the filter has its own tests.


def gene_block_overlap(*args, **kw):
    kw.setdefault("min_block_snp", 0)
    kw.setdefault("min_block_kb", 0)
    return _gene_block_overlap(*args, **kw)



def _blocks():
    # groups: A = {s1, s2}, B = {s3, s4}. Every sample appears (so all are "analyzed").
    # gene of interest: chr1 [1000, 2000].
    rows = [
        # IBD (different == 0)
        ("s1", "s2", "1", 500, 1500, 0),    # within A, overlaps gene
        ("s1", "s3", "1", 1800, 2500, 0),   # A-B, overlaps gene
        ("s2", "s4", "1", 900, 2100, 0),    # A-B, overlaps gene
        ("s1", "s2", "1", 4000, 5000, 0),   # within A, far from gene (already counted anyway)
        ("s3", "s4", "1", 6000, 7000, 0),   # within B, far from gene -> B,B numerator 0
        # a non-IBD segment (different == 1) just to prove s3/s4 are "analyzed"
        ("s3", "s4", "1", 100, 400, 1),
    ]
    return pd.DataFrame(rows, columns=["sample1", "sample2", "chr", "start", "end", "different"])


_S2G = {"s1": "A", "s2": "A", "s3": "B", "s4": "B"}
_GENES = pd.DataFrame({"name": ["g"], "chr": ["1"], "start": [1000], "end": [2000]})


def _frac(df, a, b):
    r = df[(df.group_a == a) & (df.group_b == b)].iloc[0]
    return r.n_pairs_ibd, r.n_pairs_total, r.frac_pairs_ibd


def test_overlap_counts_distinct_pairs_over_all_pairs():
    out = gene_block_overlap(_blocks(), _GENES, _S2G)
    # denominators: within A = C(2,2)=1, within B = 1, between A-B = 2*2 = 4
    assert _frac(out, "A", "A") == (1, 1, 1.0)          # s1-s2 overlaps
    assert _frac(out, "A", "B") == (2, 4, 0.5)          # s1-s3 and s2-s4 overlap (of 4 pairs)
    assert _frac(out, "B", "B") == (0, 1, 0.0)          # s3-s4 IBD only far away


def test_a_pair_with_two_overlapping_blocks_counts_once():
    b = _blocks()
    # add a SECOND s1-s2 block that also overlaps the gene
    extra = pd.DataFrame([("s1", "s2", "1", 1100, 1200, 0)],
                         columns=b.columns)
    out = gene_block_overlap(pd.concat([b, extra], ignore_index=True), _GENES, _S2G)
    assert _frac(out, "A", "A")[0] == 1                 # still one distinct pair, not two


def test_within_padding_pulls_in_a_nearby_block():
    b = pd.DataFrame([
        ("s1", "s2", "1", 500, 900, 0),      # within A, ends 100 bp before the gene
        ("s3", "s4", "1", 10, 20, 1),        # keep s3/s4 analyzed
        ("s1", "s3", "1", 10, 20, 1),        # keep the A-B pairs analyzed but non-IBD
        ("s2", "s4", "1", 10, 20, 1),
    ], columns=["sample1", "sample2", "chr", "start", "end", "different"])
    strict = gene_block_overlap(b, _GENES, _S2G, within=0)
    padded = gene_block_overlap(b, _GENES, _S2G, within=200)
    assert _frac(strict, "A", "A")[0] == 0              # 900 < 1000, no overlap
    assert _frac(padded, "A", "A")[0] == 1              # gene padded to [800, 2200]


def test_chromosome_names_are_normalised():
    b = _blocks().copy()
    b["chr"] = "Pf3D7_01_v3"                            # block uses the long spelling
    genes = pd.DataFrame({"name": ["g"], "chr": ["1"], "start": [1000], "end": [2000]})
    out = gene_block_overlap(b, genes, _S2G)
    assert _frac(out, "A", "A") == (1, 1, 1.0)          # still matches after normalise_chr


def test_gene_id_is_carried_through_when_present():
    genes = pd.DataFrame({"name": ["pfvar", "pfvar"], "gene_id": ["PF3D7_0100100", "PF3D7_0100300"],
                          "chr": ["1", "1"], "start": [1000, 5000], "end": [2000, 6000]})
    b = pd.DataFrame([
        ("s1", "s2", "1", 1500, 1800, 0),   # overlaps first pfvar only
        ("s3", "s4", "1", 10, 20, 1),
    ], columns=["sample1", "sample2", "chr", "start", "end", "different"])
    out = gene_block_overlap(b, genes, _S2G)
    assert "gene_id" in out.columns
    # the two same-named genes stay distinguishable by gene_id
    assert set(out["gene_id"].unique()) == {"PF3D7_0100100", "PF3D7_0100300"}
