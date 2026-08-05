"""Interval conventions: 0-based half-open everywhere, VCF positions 1-based."""

import numpy as np
import pandas as pd

from plasgenomicsutils.lib.intervals import (
    blocks_to_half_open,
    overlap_span,
    overlaps,
    variant_in_interval,
)
from plasgenomicsutils.lib.ibd_matrix import read_blocks


def test_blocks_to_half_open_shifts_only_end():
    df = pd.DataFrame({"chr": ["c1"], "start": [100], "end": [200]})
    out = blocks_to_half_open(df)
    assert out["start"].tolist() == [100]
    assert out["end"].tolist() == [201]          # inclusive last SNP -> exclusive end
    assert df["end"].tolist() == [200]           # input untouched


def test_read_blocks_returns_half_open(tmp_path):
    p = tmp_path / "b.tsv"
    p.write_text(
        "sample1\tsample2\tchr\tstart\tend\tdifferent\tNsnp\n"
        "a\tb\tPf3D7_01_v3\t100\t200\t0\t5\n"
    )
    df = read_blocks(str(p))
    assert df["start"].iloc[0] == 100
    assert df["end"].iloc[0] == 201
    # the segment spans the last SNP, so its width counts that base
    assert df["end"].iloc[0] - df["start"].iloc[0] == 101


def test_overlaps_is_half_open():
    starts = np.array([0, 100, 200, 300])
    ends = np.array([50, 200, 300, 400])
    # [100,200) vs [200,300): abutting intervals do NOT overlap
    m = overlaps(starts, ends, 200, 300)
    assert m.tolist() == [False, False, True, False]
    # a one-base query picks up only the interval containing it
    assert overlaps(starts, ends, 150, 151).tolist() == [False, True, False, False]


def test_overlap_span_clips_to_intersection():
    assert overlap_span(100, 300, 200, 400) == (200, 300)
    assert overlap_span(100, 300, 100, 300) == (100, 300)
    lo, hi = overlap_span(100, 150, 200, 300)     # disjoint -> empty (lo >= hi)
    assert lo >= hi


def test_variant_in_interval_is_zero_based():
    # a variant sits in [100, 200) when its 0-based position does; the end is exclusive
    pos = np.array([99, 100, 199, 200])
    assert variant_in_interval(pos, 100, 200).tolist() == [False, True, True, False]


def test_snp_label_is_the_only_id_builder():
    from plasgenomicsutils.lib.intervals import snp_label, vcf_pos
    assert snp_label("Pf3D7_07_v3", 403221) == "Pf3D7_07_v3:403221"
    assert snp_label(["c1", "c2"], [5, 6]) == ["c1:5", "c2:6"]
    assert vcf_pos(403221) == 403222          # 1-based only for human cross-reference


def test_unstamped_af_table_is_rejected(tmp_path):
    """The silent-corruption case: a dense panel shifted by one still *matches*.

    299 of 300 labels find a real neighbouring SNP, so only one looks 'missing' and the
    join would otherwise succeed with every AF attached to the wrong SNP. The stamp is
    what makes this detectable.
    """
    import gzip
    import pytest
    from plasgenomicsutils.lib.ibd_selection import load_global_af

    labels = [f"chr1:{p}" for p in range(100, 400)]
    old = pd.DataFrame({"snp_id": [f"chr1:{p + 1}" for p in range(100, 400)], "af": 0.3})
    path = tmp_path / "af.tsv.gz"
    with gzip.open(path, "wt") as fh:
        old.to_csv(fh, sep="\t", index=False)
    with pytest.raises(SystemExit, match="coordinate-system stamp"):
        load_global_af(str(path), labels)


def test_stamped_af_table_round_trips(tmp_path):
    from plasgenomicsutils.lib.intervals import SNP_COORD_SYSTEM
    from plasgenomicsutils.lib.ibd_selection import load_global_af
    from plasgenomicsutils.utils.small_utils import Utils

    labels = ["chr1:99", "chr1:199"]
    af = pd.DataFrame({"snp_id": labels, "af": [0.25, 0.75]})
    path = tmp_path / "af.tsv.gz"
    Utils.write_tsv_gz(af, str(path), header_comment=f"snp_coord_system={SNP_COORD_SYSTEM}")
    assert load_global_af(str(path), labels).tolist() == [0.25, 0.75]


def test_source_ids_are_kept_but_not_used_as_keys(tmp_path):
    """A BED name column written with %POS (1-based) must not become the key."""
    from plasgenomicsutils.lib.vcf_io import SnpPanel

    bed = tmp_path / "s.bed"
    bed.write_text("chr1\t99\t100\tchr1:100\nchr1\t199\t200\tchr1:200\n")
    p = SnpPanel.from_bed(str(bed)).df
    assert list(p.snp_id) == ["chr1:99", "chr1:199"]      # derived from the 0-based start
    assert list(p.source_id) == ["chr1:100", "chr1:200"]  # the file's own id, kept aside
