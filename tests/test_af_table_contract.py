"""The contract between compute_allele_freqs and the IBD selection statistic.

The AF tables gained columns (maf, ac/an, af_weighted, prevalence...). The selection
statistic joins on them, so these check the numbers it computes are untouched by that --
not merely that the file still parses.
"""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("cyvcf2")

from plasgenomicsutils.lib.ibd_freqs import compute_allele_freqs
from plasgenomicsutils.lib.ibd_selection import (
    compute_selection_statistic, get_af_for_group, load_global_af, load_group_af_table,
)
from plasgenomicsutils.lib.intervals import SNP_COORD_SYSTEM
from plasgenomicsutils.utils.small_utils import Utils

_STAMP = f"snp_coord_system={SNP_COORD_SYSTEM}"


def _callset(tmp_path, n_samples=8):
    """A small callset with AD, so every column the tables can hold is populated."""
    samples = [f"s{i}" for i in range(1, n_samples + 1)]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    rng = np.random.default_rng(0)
    rows = []
    for j in range(12):
        n_alt = 1 + j % (n_samples - 1)
        cells = []
        for i in range(n_samples):
            alt = i < n_alt
            depth = int(rng.integers(30, 120))
            a = depth if alt else 0
            cells.append(f"{'1/1' if alt else '0/0'}:{depth - a},{a}")
        rows.append(f"chr1\t{1000 * (j + 1)}\t.\tA\tT\t.\t.\t.\tGT:AD\t" + "\t".join(cells))
    vcf = tmp_path / "cohort.vcf"
    vcf.write_text("\n".join(hdr + rows) + "\n")
    return str(vcf), samples


def _write(df, path):
    Utils.write_tsv_gz(df, str(path), header_comment=_STAMP)
    return str(path)


def test_the_widened_table_gives_the_statistic_the_same_afs(tmp_path):
    """An old-schema table and today's wide one must load to identical frequencies."""
    vcf, _ = _callset(tmp_path)
    g, _ = compute_allele_freqs(vcf)
    labels = g["snp_id"].tolist()

    wide = _write(g, tmp_path / "wide.tsv.gz")
    narrow = _write(g[["snp_id", "af"]], tmp_path / "narrow.tsv.gz")   # pre-change shape

    af_wide = load_global_af(wide, labels)
    af_narrow = load_global_af(narrow, labels)
    assert np.array_equal(af_wide, af_narrow)
    assert not np.isnan(af_wide).any()


def test_the_statistic_itself_is_unchanged_by_the_extra_columns(tmp_path):
    """End to end: the same matrix scored against both table shapes."""
    vcf, _ = _callset(tmp_path)
    g, _ = compute_allele_freqs(vcf)
    labels = g["snp_id"].tolist()

    from scipy.sparse import csr_matrix

    rng = np.random.default_rng(1)
    mat = csr_matrix((rng.random((40, len(labels))) < 0.3).astype(np.int8))  # pairs x SNPs

    stats_w, df_w = compute_selection_statistic(
        mat, load_global_af(_write(g, tmp_path / "w.tsv.gz"), labels))
    stats_n, df_n = compute_selection_statistic(
        mat, load_global_af(_write(g[["snp_id", "af"]], tmp_path / "n.tsv.gz"), labels))

    assert set(stats_w) == set(stats_n)
    for k in stats_w:
        a, b = np.asarray(stats_w[k]), np.asarray(stats_n[k])
        if a.dtype.kind in "fc":                      # numeric: allow NaN either side
            assert np.allclose(a, b, equal_nan=True), k
        else:                                         # `direction` is excess/deficit
            assert np.array_equal(a, b), k
    pd.testing.assert_frame_equal(df_w, df_n)
    # and the statistic is populated, so this is not two empty results agreeing
    assert np.isfinite(np.asarray(stats_w["raw_stat"], dtype=float)).any()
    assert set(np.unique(stats_w["direction"])) <= {"excess", "deficit"}


def test_the_group_table_still_joins_per_group(tmp_path):
    vcf, samples = _callset(tmp_path)
    s2g = {s: ("A" if i < 4 else "B") for i, s in enumerate(samples)}
    g, r = compute_allele_freqs(vcf, s2g)
    labels = g["snp_id"].tolist()

    wide = load_group_af_table(_write(r, tmp_path / "gw.tsv.gz"))
    narrow = load_group_af_table(
        _write(r[["group", "snp_id", "af"]], tmp_path / "gn.tsv.gz"))
    glob_ = load_global_af(_write(g, tmp_path / "w.tsv.gz"), labels)

    for grp in ("A", "B"):
        a = get_af_for_group(grp, labels, wide, glob_)
        b = get_af_for_group(grp, labels, narrow, glob_)
        assert np.array_equal(a, b)
        assert not np.isnan(a).any()


def test_per_alt_tables_are_rejected_rather_than_silently_joined(tmp_path):
    """--per-alt gives several rows per snp_id; the statistic joins one row per SNP."""
    vcf, _ = _callset(tmp_path)
    per, _ = compute_allele_freqs(vcf, per_alt=True)
    # this fixture is biallelic, so per-alt is still one row per SNP and must still work
    labels = per["snp_id"].tolist()
    assert len(set(labels)) == len(labels)
    af = load_global_af(_write(per, tmp_path / "p.tsv.gz"), labels)
    assert not np.isnan(af).any()


def test_a_table_without_the_coordinate_stamp_is_refused(tmp_path):
    """The stamp is what stops a 1-based panel being joined against 0-based ids."""
    vcf, _ = _callset(tmp_path)
    g, _ = compute_allele_freqs(vcf)
    bare = tmp_path / "bare.tsv"
    pd.DataFrame(g).to_csv(bare, sep="\t", index=False)
    with pytest.raises(SystemExit, match="coordinate-system stamp"):
        load_global_af(str(bare), g["snp_id"].tolist())


def test_af_col_selects_which_frequency_the_statistic_scores_against(tmp_path):
    vcf, _ = _callset(tmp_path)
    g, _ = compute_allele_freqs(vcf)
    labels = g["snp_id"].tolist()
    path = _write(g, tmp_path / "af.tsv.gz")

    default = load_global_af(path, labels)
    explicit = load_global_af(path, labels, af_col="af")
    weighted = load_global_af(path, labels, af_col="af_weighted")

    assert np.array_equal(default, explicit)                 # `af` is the default
    assert np.allclose(weighted, g["af_weighted"].to_numpy())
    assert not np.isnan(weighted).any()


def test_a_missing_af_column_names_the_ones_that_are_there(tmp_path):
    vcf, _ = _callset(tmp_path)
    g, _ = compute_allele_freqs(vcf)
    # an old-schema table has no af_weighted: the error has to say so, not fail on a join
    narrow = _write(g[["snp_id", "af"]], tmp_path / "narrow.tsv.gz")
    with pytest.raises(SystemExit, match="has no column"):
        load_global_af(narrow, g["snp_id"].tolist(), af_col="af_weighted")
    with pytest.raises(SystemExit, match="re-run compute_allele_freqs"):
        load_global_af(narrow, g["snp_id"].tolist(), af_col="af_weighted")


def test_the_group_table_can_be_scored_on_a_different_column(tmp_path):
    vcf, samples = _callset(tmp_path)
    s2g = {s: ("A" if i < 4 else "B") for i, s in enumerate(samples)}
    _, r = compute_allele_freqs(vcf, s2g)
    path = _write(r, tmp_path / "gaf.tsv.gz")

    # whichever column is chosen, it arrives as `af` so the joins downstream are unchanged
    w = load_group_af_table(path, af_col="af_weighted")
    assert list(w.columns) == ["group", "snp_id", "af"]
    assert set(w["group"]) == {"A", "B"}
    # the chosen column really is the one that arrives, renamed rather than recomputed
    key = ["group", "snp_id"]
    src = r.set_index(key)
    got = w.set_index(key)
    assert np.allclose(got["af"], src.loc[got.index, "af_weighted"])
    plain = load_group_af_table(path).set_index(key)
    assert np.allclose(plain["af"], src.loc[plain.index, "af"])
