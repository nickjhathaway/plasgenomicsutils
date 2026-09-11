"""LD decay. The heavy pairwise work is scikit-allel's rogers_huff_r; these tests pin the
windowing, binning and grouping around it, and that missing calls cost a pair rather than
a SNP."""

from __future__ import annotations

import numpy as np
import pytest

from plasgenomicsutils.lib.ld import _pairs_one_chrom, _thin, half_decay, ld_decay


def _linked(n_samples=60, n_var=12, noise=0.0, seed=0):
    """Variants along a chromosome, each copying its neighbour with a per-step flip
    probability. Correlation then decays as (1 - 2*noise)^steps, so r-squared genuinely
    falls with distance -- a flip probability that grows past 0.5 instead produces
    *anti*-correlation and r-squared climbing back up, which is not decay."""
    rng = np.random.default_rng(seed)
    cur = rng.integers(0, 2, n_samples)
    rows = [cur.copy()]
    for _ in range(n_var - 1):
        cur = np.where(rng.random(n_samples) < noise, 1 - cur, cur)
        rows.append(cur.copy())
    return np.array(rows, dtype=np.int8) * 2


def test_perfectly_linked_variants_have_r2_of_one():
    gn = _linked(noise=0.0)
    pos = np.arange(gn.shape[0]) * 1000
    d, r2 = _pairs_one_chrom(gn, pos, max_dist=50_000)
    assert len(d) == gn.shape[0] * (gn.shape[0] - 1) // 2      # every pair is in range
    assert r2 == pytest.approx(1.0)


def test_only_pairs_within_max_dist_are_counted():
    gn = _linked(n_var=6)
    pos = np.array([0, 1000, 2000, 50_000, 51_000, 52_000])
    d, _ = _pairs_one_chrom(gn, pos, max_dist=5000)
    # three pairs inside each cluster, none across the 48 kb gap
    assert len(d) == 6
    assert d.max() <= 5000


def test_blocking_does_not_change_the_pair_set(monkeypatch):
    """Blocks exist only to bound memory; each pair must be counted exactly once."""
    import plasgenomicsutils.lib.ld as ld

    gn = _linked(n_var=40, noise=0.01)
    pos = np.arange(gn.shape[0]) * 1000
    d_ref, r_ref = _pairs_one_chrom(gn, pos, max_dist=10_000)
    monkeypatch.setattr(ld, "_BLOCK", 5)                        # force many small blocks
    d_small, r_small = ld._pairs_one_chrom(gn, pos, max_dist=10_000)
    order_ref, order_small = np.argsort(d_ref, kind="stable"), np.argsort(d_small, kind="stable")
    assert len(d_ref) == len(d_small)
    assert np.array_equal(d_ref[order_ref], d_small[order_small])
    assert r_ref[order_ref] == pytest.approx(r_small[order_small], abs=1e-6)


def test_a_missing_call_costs_the_pair_not_the_snp():
    gn = _linked(n_var=4)
    pos = np.arange(4) * 1000
    full = _pairs_one_chrom(gn, pos, max_dist=50_000)[1]
    gm = gn.copy()
    gm[0, :20] = -1                                             # a third of one SNP
    part = _pairs_one_chrom(gm, pos, max_dist=50_000)[1]
    assert len(part) == len(full)                               # the SNP is still scanned
    assert np.all(np.isfinite(part))


def test_decay_falls_with_distance_and_bins_are_half_open():
    gn = _linked(n_var=40, noise=0.05, seed=3)
    pos = np.arange(gn.shape[0]) * 1000
    chrom = np.array(["c1"] * gn.shape[0])
    df, half = ld_decay(gn, chrom, pos, max_dist=40_000, bins=4, maf=0.0)
    ok = df.dropna(subset=["mean_r2"])
    assert ok["mean_r2"].iloc[0] > ok["mean_r2"].iloc[-1]
    # bins tile [0, max_dist) with no gaps and no overlap
    assert list(ok["bin_start"]) == [0, 10_000, 20_000, 30_000]
    assert list(ok["bin_end"]) == [10_000, 20_000, 30_000, 40_000]
    assert int(df["n_pairs"].sum()) == len(_pairs_one_chrom(gn, pos, 40_000)[0])


def test_groups_are_scanned_separately_and_small_ones_skipped():
    gn = np.hstack([_linked(n_samples=40, n_var=10, noise=0.0),          # group a: linked
                    _linked(n_samples=40, n_var=10, noise=0.5, seed=7),  # group b: not
                    _linked(n_samples=2, n_var=10)])                     # too small
    pos = np.arange(10) * 1000
    chrom = np.array(["c1"] * 10)
    groups = np.array(["a"] * 40 + ["b"] * 40 + ["tiny"] * 2)
    df, _ = ld_decay(gn, chrom, pos, groups=groups, max_dist=20_000, bins=2, maf=0.0)
    assert set(df["group"]) == {"a", "b"}                       # 'tiny' is below min_samples
    a = df[df.group == "a"]["mean_r2"].iloc[0]
    b = df[df.group == "b"]["mean_r2"].iloc[0]
    assert a > b


def test_thinning_is_even_and_bounded():
    idx = np.arange(1000)
    assert np.array_equal(_thin(idx, 5000), idx)                # under budget: untouched
    t = _thin(idx, 100)
    assert len(t) == 100
    assert t[0] == 0 and t[-1] == 999                           # spans the chromosome
    assert np.all(np.diff(t) > 0)


def test_half_decay_interpolates_and_reports_nothing_when_flat():
    import pandas as pd

    df = pd.DataFrame({"group": ["a"] * 3 + ["b"] * 3,
                       "bin_mid": [1000, 2000, 3000] * 2,
                       "mean_r2": [0.4, 0.3, 0.1, 0.2, 0.2, 0.2]})
    h = half_decay(df).set_index("group")["half_decay_bp"]
    assert 2000 < h["a"] < 3000        # crosses 0.2 between the second and third bins
    assert np.isnan(h["b"])            # never halves


# --- multiallelic sites: decision 2 is to drop them, counted -------------------------

def _ld_vcf(tmp_path):
    """Two biallelic SNPs and one triallelic, all with the same carriers.

    `gt_types` reports 3 for any homozygous-alternate call, so `1/1` and `2/2` were both
    coded as dosage 2 -- two different alleles collapsed into one symbol. That inflates r²
    between multiallelic sites, and unlike a dropped site it leaves no trace in the SNP
    count.
    """
    samples = [f"s{i}" for i in range(1, 9)]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    rows = [
        ("1000", "C", ["0/0"] * 4 + ["1/1"] * 4),
        ("2000", "C,G", ["0/0"] * 4 + ["1/1"] * 2 + ["2/2"] * 2),
        ("3000", "C", ["0/0"] * 4 + ["1/1"] * 4),
    ]
    for pos, alt, gts in rows:
        hdr.append(f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts))
    p = tmp_path / "ld.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def test_multiallelic_sites_are_dropped_not_merged(tmp_path):
    from plasgenomicsutils.lib.ld import read_dosages

    gn, chrom, pos, _names, st = read_dosages(_ld_vcf(tmp_path))
    assert gn.shape[0] == 2, "the triallelic record must not reach the dosage matrix"
    assert list(pos) == [999, 2999]
    assert st["multiallelic_skipped"] == 1
    assert st["variants_read"] == 3


def test_the_drop_is_reported_rather_than_silent(tmp_path, capsys):
    """A dropped site is visible in the SNP count; a merged one is not."""
    from plasgenomicsutils.lib.ld import read_dosages

    read_dosages(_ld_vcf(tmp_path))
    out = capsys.readouterr().out
    assert "multiallelic" in out
    assert "1" in out


def test_a_biallelic_only_file_reports_no_skips_and_is_unchanged(tmp_path):
    from plasgenomicsutils.lib.ld import read_dosages

    samples = [f"s{i}" for i in range(1, 9)]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for pos in (1000, 2000):
        hdr.append(f"chr1\t{pos}\t.\tA\tC\t.\t.\t.\tGT\t"
                   + "\t".join(["0/0"] * 4 + ["1/1"] * 4))
    p = tmp_path / "bi.vcf"
    p.write_text("\n".join(hdr) + "\n")
    gn, _c, pos, _n, st = read_dosages(str(p))
    assert gn.shape[0] == 2
    assert st["multiallelic_skipped"] == 0
    assert list(pos) == [999, 1999]


def test_the_two_maf_definitions_now_agree(tmp_path):
    """`ld_decay --maf` and `maf_filter --maf-min` must mean the same thing.

    They did not: LD computed the *sum of alternates*, which `maf_filter`'s docstring
    explicitly rejects, while `maf_filter` uses the second-most-common allele. The two
    differ only at a multiallelic site -- so dropping those, which decision 2 asks for
    anyway, is what makes the two definitions coincide.
    """
    from plasgenomicsutils.lib.ld import _maf, read_dosages

    gn, _c, _p, _n, _st = read_dosages(_ld_vcf(tmp_path))
    maf, n = _maf(gn)
    # every surviving site is biallelic 4/4, so both readings give 0.5
    assert list(n) == [8, 8]
    assert maf.tolist() == pytest.approx([0.5, 0.5])


# --- `*` is not an alternate base -----------------------------------------------------
#
# `read_dosages` skips multiallelic records, which is right: r-squared is a squared
# correlation between two binary indicators and a multiallelic locus has no unique scalar
# summary. But it counted `*` towards that, so `A > T,*` -- a biallelic SNP with a
# spanning-deletion note attached -- was thrown away as multiallelic. On a real Uganda
# callset that removed a large share of the panel, and LD decay is precisely the analysis
# that needs SNP density.
#
# The calls that ARE the deletion still have to go: there is no base there to correlate, so
# they read as missing, which is what -1 already means to rogers_huff_r.

_LD_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
)


def _star_vcf(tmp_path, samples, rows, name="ld.vcf"):
    hdr = _LD_HDR + ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                     + "\t".join(samples) + "\n")
    body = "".join(f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts) + "\n"
                   for pos, alt, gts in rows)
    p = tmp_path / name
    p.write_text(hdr + body)
    return str(p)


def test_a_snp_beside_a_spanning_deletion_is_read_not_skipped(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.ld import read_dosages

    samples = [f"s{i}" for i in range(1, 7)]
    rows = [
        (1000, "T", ["0/0", "0/0", "0/0", "1/1", "1/1", "1/1"]),
        (2000, "T,*", ["0/0", "0/0", "2/2", "1/1", "1/1", "1/1"]),   # SNP + deletion
        (3000, "T,G", ["0/0", "1/1", "1/1", "2/2", "2/2", "0/0"]),   # truly multiallelic
    ]
    gn, _chrom, pos, _names, counts = read_dosages(_star_vcf(tmp_path, samples, rows))
    assert list(pos) == [999, 1999]
    assert counts["multiallelic_skipped"] == 1        # only the real one
    assert counts["spanning_del_masked"] == 1


def test_the_deleted_haplotype_reads_as_missing_not_reference(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.ld import read_dosages

    samples = [f"s{i}" for i in range(1, 7)]
    rows = [(2000, "T,*", ["0/0", "0/0", "2/2", "1/1", "1/1", "1/1"])]
    gn, _c, _p, _n, _k = read_dosages(_star_vcf(tmp_path, samples, rows))
    # s3 carries the deletion: -1, not 0 (which would say "confidently reference") and not
    # 2 (which would say "carries T")
    assert list(gn[0]) == [0, 0, -1, 2, 2, 2]


def test_a_record_with_nothing_but_a_deletion_carries_no_snp(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.ld import read_dosages

    samples = [f"s{i}" for i in range(1, 7)]
    rows = [
        (1000, "T", ["0/0", "0/0", "0/0", "1/1", "1/1", "1/1"]),
        (2000, "*", ["0/0", "0/0", "1/1", "1/1", "0/0", "0/0"]),
    ]
    _gn, _c, pos, _n, counts = read_dosages(_star_vcf(tmp_path, samples, rows))
    assert list(pos) == [999]
    assert counts["multiallelic_skipped"] == 1
