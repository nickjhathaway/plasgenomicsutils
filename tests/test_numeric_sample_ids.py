"""Numeric sample ids join like the labels they are.

A cohort identified by numbers -- micronix ids, barcodes, anything without a letter --
reads back from a metadata table as int64, while every id it has to be matched against (a
VCF sample name, an IBD pair label) is a string. The two never compare equal and nothing
raises: the merge simply finds nothing, so the run completes, writes every output, and
reports no groups.
"""

import pandas as pd
import pytest

from plasgenomicsutils.lib.ibd_analyze import annotate_pairs_with_groups
from plasgenomicsutils.utils.small_utils import Utils


def _meta_file(tmp_path, ids, name="meta.tsv", region=("East", "North", "West")):
    rows = ["sample\tregion"] + [f"{s}\t{region[i % len(region)]}" for i, s in enumerate(ids)]
    p = tmp_path / name
    p.write_text("\n".join(rows) + "\n")
    return str(p)


def _pairs(ids):
    return pd.DataFrame([{"sample1": str(a), "sample2": str(b), "frac_ibd": 0.1}
                         for i, a in enumerate(ids) for b in ids[i + 1:]])


NUMERIC = [4064928010, 4074600944, 4074608135]


def test_all_numeric_sample_ids_are_read_as_strings(tmp_path):
    meta = Utils.read_meta(_meta_file(tmp_path, NUMERIC), wants=("sample", "region"),
                           quiet=True)
    assert pd.api.types.is_string_dtype(meta["sample"])
    assert meta["sample"].tolist() == [str(x) for x in NUMERIC]


def test_a_numeric_cohort_still_forms_groups(tmp_path):
    """The regression this file exists for: every pair came back 'unknown'."""
    meta = Utils.read_meta(_meta_file(tmp_path, NUMERIC), wants=("sample", "region"),
                           quiet=True)
    ann = annotate_pairs_with_groups(_pairs(NUMERIC), meta, "region")
    assert "unknown" not in set(ann["group1"]) | set(ann["group2"])
    assert sorted(set(ann["group1"]) | set(ann["group2"])) == ["East", "North", "West"]


def test_a_blank_cell_does_not_turn_ids_into_floats(tmp_path):
    """One missing value makes pandas read the column as float64, and `4064928010.0`
    matches nothing either."""
    p = tmp_path / "gappy.tsv"
    p.write_text("sample\tregion\n4064928010\tEast\n\tNorth\n4074608135\tWest\n")
    meta = Utils.read_meta(str(p), wants=("sample", "region"), quiet=True)
    assert meta["sample"].tolist() == ["4064928010", "", "4074608135"]


def test_string_ids_are_left_exactly_as_they_are(tmp_path):
    ids = ["RCN13010", "RCN13048"]
    meta = Utils.read_meta(_meta_file(tmp_path, ids), wants=("sample", "region"), quiet=True)
    assert meta["sample"].tolist() == ids


def test_annotate_coerces_even_a_hand_built_frame(tmp_path):
    """The function is called directly too, so it cannot rely on read_meta having run."""
    meta = pd.DataFrame({"sample": NUMERIC, "region": ["East", "North", "West"]})
    ann = annotate_pairs_with_groups(_pairs(NUMERIC), meta, "region")
    assert "unknown" not in set(ann["group1"]) | set(ann["group2"])


def test_metadata_that_matches_nothing_says_so_loudly(capsys):
    """After the coercion, no overlap means a genuinely mismatched file -- which must not
    look like a clean run that happened to find no groups."""
    meta = pd.DataFrame({"sample": ["OTHER-1", "OTHER-2"], "region": ["East", "North"]})
    ann = annotate_pairs_with_groups(_pairs(NUMERIC), meta, "region")
    out = capsys.readouterr().out
    assert "no sample in the metadata matches" in out
    assert "4064928010" in out and "OTHER-1" in out          # examples from both sides
    assert set(ann["group1"]) == {"unknown"}


def test_a_partial_match_is_not_warned_about(capsys):
    """Some samples missing from the metadata is ordinary; only zero overlap is a fault."""
    meta = pd.DataFrame({"sample": [str(NUMERIC[0])], "region": ["East"]})
    annotate_pairs_with_groups(_pairs(NUMERIC), meta, "region")
    assert "no sample in the metadata matches" not in capsys.readouterr().out
