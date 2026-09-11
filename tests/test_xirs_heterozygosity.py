"""The IBD selection statistic's two uses of allele frequency, generalised to k alleles.

Allele frequency enters XiR,s in exactly two places: a scale divisor `sqrt(p(1-p))` and the
binning variable. Both are heterozygosity terms wearing a biallelic disguise, since
`p(1-p) == H/2` where `H = 1 - sum(p_i^2)`. So both generalise by substitution, and on a
biallelic panel the substitution changes nothing at all -- which is the property these tests
exist to hold onto.

The statistic itself never compares alleles. It counts IBD-interval coverage per SNP, and the
allele-comparison problem lives upstream in hmmibd's emission model. So there was never
anything here that needed replacing.
"""

import numpy as np
import pytest
from scipy import sparse

from plasgenomicsutils.lib.ibd_selection import compute_selection_statistic


def _matrix(n_pairs=400, n_snps=300, seed=0):
    rng = np.random.default_rng(seed)
    m = (rng.random((n_pairs, n_snps)) < 0.25).astype(np.int8)
    return sparse.csr_matrix(m)


def _biallelic_he(af):
    return 2.0 * af * (1.0 - af)


def test_on_a_biallelic_panel_the_heterozygosity_form_changes_nothing(tmp_path):
    """`sqrt(he/2)` is `sqrt(p(1-p))` and binning on `he` is binning on `maf`.

    Both are exact, not approximate: `he = 2p(1-p)` at k=2, and `he` is strictly monotone in
    `maf` on [0, 0.5], so equal-frequency bins over the two are the same bins. Any drift here
    means an existing biallelic result moved, which is the one thing this must not do.
    """
    rng = np.random.default_rng(1)
    af = rng.uniform(0.05, 0.95, 300)
    mat = _matrix(n_snps=300)

    old, old_bins = compute_selection_statistic(mat, af, n_bins=20)
    new, new_bins = compute_selection_statistic(mat, af, n_bins=20, he=_biallelic_he(af))

    np.testing.assert_allclose(old["raw_stat"], new["raw_stat"], rtol=0, atol=0)
    np.testing.assert_allclose(old["z_score"], new["z_score"], rtol=0, atol=0)
    assert list(old_bins["n_snps"]) == list(new_bins["n_snps"])


def test_the_scale_term_is_the_heterozygosity_one():
    """`sqrt(he/2)` at k=2 *is* `sqrt(p(1-p))`, and at k=3 it is larger than `p(1-p)` says."""
    from plasgenomicsutils.lib.ibd_selection import _scale_and_gate

    af = np.array([0.2, 0.5, 0.8])
    _v, biallelic = _scale_and_gate(af, None)
    _v, from_he = _scale_and_gate(af, _biallelic_he(af))
    np.testing.assert_allclose(biallelic, from_he, rtol=0, atol=0)

    # 0.50/0.25/0.25 -> he = 0.625, so the scale is sqrt(0.3125) where the collapsed
    # af = 0.5 would have said sqrt(0.25): a 12% understatement of the divisor
    _v, tri = _scale_and_gate(np.array([0.5]), np.array([0.625]))
    assert tri[0] == pytest.approx(np.sqrt(0.3125))
    assert tri[0] > np.sqrt(0.25)


def test_a_reference_absent_site_is_kept_rather_than_silently_deleted():
    """`af == 1` means every sample carries an alternate, not that the site is monomorphic.

    The old gate was `0 < af < 1`, which deleted exactly the maximally informative sites: a
    4C/4G site is perfectly polymorphic (`he = 0.5`) and collapses to `af = 1.0`. Nothing said
    so, and the Bonferroni denominator shrank with it.
    """
    af = np.array([0.4, 1.0, 0.6, 0.0, 0.5])
    he = np.array([0.48, 0.5, 0.48, 0.5, 0.5])      # both the af==1 and af==0 sites are real
    mat = _matrix(n_pairs=60, n_snps=5)

    old, _ = compute_selection_statistic(mat, af, n_bins=1)
    new, _ = compute_selection_statistic(mat, af, n_bins=1, he=he)
    assert np.isnan(old["raw_stat"][1]) and np.isnan(old["raw_stat"][3])
    assert np.isfinite(new["raw_stat"][1]) and np.isfinite(new["raw_stat"][3])


def test_a_truly_monomorphic_site_is_still_dropped():
    """`he == 0` is the honest test for "nothing to see here"."""
    af = np.array([0.4, 0.5, 0.6])
    he = np.array([0.48, 0.0, 0.48])
    mat = _matrix(n_pairs=60, n_snps=3)
    st, _ = compute_selection_statistic(mat, af, n_bins=1, he=he)
    assert np.isnan(st["raw_stat"][1])
    assert np.isfinite(st["raw_stat"][0]) and np.isfinite(st["raw_stat"][2])


def test_binning_puts_a_multiallelic_site_above_every_biallelic_one():
    """A triallelic site can be more informative than any biallelic one, and must bin there.

    `maf` caps at 0.5; `1 - max(p_i)` reaches `1 - 1/k`. Binning a 0.5/0.25/0.25 site with the
    most balanced biallelic ones compares it against a *less* informative reference class.

    Note this is `maf_k`, not `he`: heterozygosity would order them the same way but is flat
    near 0.5, which moves biallelic sites between bins -- see
    `test_the_binning_key_is_unchanged_on_a_tied_biallelic_panel`.
    """
    af = np.concatenate([np.linspace(0.05, 0.95, 40), [0.5]])
    n_alleles = np.concatenate([np.full(40, 2), [3]])
    maf_k = np.concatenate([np.where(np.linspace(0.05, 0.95, 40) <= 0.5,
                                     np.linspace(0.05, 0.95, 40),
                                     1 - np.linspace(0.05, 0.95, 40)), [2 / 3]])
    mat = _matrix(n_pairs=80, n_snps=41)
    _st, bins = compute_selection_statistic(mat, af, n_bins=4,
                                            n_alleles=n_alleles, maf_k=maf_k)
    assert bins["binkey_max"].iloc[-1] == pytest.approx(2 / 3, abs=1e-6)


def test_the_bin_table_says_which_variable_it_binned_on():
    """A reader has to be able to tell a k-allele-binned run from a plain MAF one."""
    af = np.linspace(0.1, 0.9, 30)
    mat = _matrix(n_pairs=60, n_snps=30)
    _s, a = compute_selection_statistic(mat, af, n_bins=3)
    _s, b = compute_selection_statistic(mat, af, n_bins=3, n_alleles=np.full(30, 2),
                                        maf_k=np.minimum(af, 1 - af))
    assert "maf_min" in a.columns and "binkey_min" not in a.columns
    assert "binkey_min" in b.columns and "maf_min" not in b.columns


def test_the_chunked_path_agrees_with_the_dense_one():
    """The two implementations must not drift, and `he` is threaded through both."""
    from plasgenomicsutils.lib.ibd_selection import _compute_chunked, _compute_dense

    rng = np.random.default_rng(3)
    af = rng.uniform(0.05, 0.95, 60)
    he = _biallelic_he(af)
    mat = _matrix(n_pairs=100, n_snps=60, seed=4)
    n_alleles = np.full(len(af), 2)
    d, _ = _compute_dense(mat, af, 5, he=he, n_alleles=n_alleles)
    c, _ = _compute_chunked(mat, af, 5, chunk_size=17, he=he, n_alleles=n_alleles)
    np.testing.assert_allclose(d["raw_stat"], c["raw_stat"], rtol=1e-12)


def _af_table(tmp_path, with_he=True, n=40):
    """A small AF table shaped like `compute_allele_freqs` writes."""
    import pandas as pd

    from plasgenomicsutils.lib.intervals import SNP_COORD_SYSTEM
    from plasgenomicsutils.utils.small_utils import Utils

    rng = np.random.default_rng(7)
    af = rng.uniform(0.05, 0.95, n)
    df = pd.DataFrame({"snp_id": [f"chr1:{1000 * (i + 1)}" for i in range(n)], "af": af})
    if with_he:
        df["he"] = _biallelic_he(af)
    p = tmp_path / ("with_he.tsv.gz" if with_he else "no_he.tsv.gz")
    Utils.write_tsv_gz(df, str(p), header_comment=f"snp_coord_system={SNP_COORD_SYSTEM}")
    return str(p), df["snp_id"].tolist()


def test_he_is_read_when_the_table_has_it(tmp_path):
    from plasgenomicsutils.lib.ibd_selection import load_global_he

    path, labels = _af_table(tmp_path)
    he = load_global_he(path, labels)
    assert he is not None and len(he) == len(labels)
    assert np.isfinite(he).all()


def test_a_table_written_before_the_column_existed_still_works(tmp_path):
    """`None` rather than an error: the biallelic forms are what such a table requires."""
    from plasgenomicsutils.lib.ibd_selection import load_global_he

    path, labels = _af_table(tmp_path, with_he=False)
    assert load_global_he(path, labels) is None


def test_a_per_alt_table_is_refused_here_too(tmp_path):
    """The `he` loader joins by snp_id like the others, so it needs the same guard."""
    import pandas as pd

    from plasgenomicsutils.lib.ibd_selection import load_global_he
    from plasgenomicsutils.lib.intervals import SNP_COORD_SYSTEM
    from plasgenomicsutils.utils.small_utils import Utils

    df = pd.DataFrame({"snp_id": ["chr1:1000", "chr1:1000", "chr1:2000"],
                       "he": [0.5, 0.5, 0.4]})
    p = tmp_path / "per_alt.tsv.gz"
    Utils.write_tsv_gz(df, str(p), header_comment=f"snp_coord_system={SNP_COORD_SYSTEM}")
    with pytest.raises(SystemExit, match="per-alt"):
        load_global_he(str(p), ["chr1:1000", "chr1:2000"])


# --- the binning key, and why it is not `he` ----------------------------------------

def _tied_af(n=4000, an=1036, seed=11):
    """Allele frequencies as a real callset has them: counts over a fixed allele number.

    The Phase 5 test used `rng.uniform`, which produces 4000 distinct values and no ties. A
    real panel has `ac/an` with `an` around a thousand, so thousands of SNPs share each
    frequency -- and ties are exactly what the binning is sensitive to.
    """
    rng = np.random.default_rng(seed)
    ac = rng.integers(1, an, n)
    return ac / an


def test_the_binning_key_is_unchanged_on_a_tied_biallelic_panel():
    """Real allele frequencies tie heavily, and `he` is flat near 0.5.

    Two distinct MAFs either side of 0.5 -- 0.48 and 0.52, say -- have *the same*
    heterozygosity. So binning on `he` merges frequency classes that binning on `maf` keeps
    apart, and an equal-frequency cut then falls in a different place. On the real 27k-SNP
    Uganda panel that moved 248 SNPs between bins, and because a bin's mean and sd shift with
    its membership, **12,986 z-scores changed**, by up to 4.6.

    So the binning key must reduce to `maf` exactly at a biallelic site, which `he` does not.
    """
    from plasgenomicsutils.lib.ibd_selection import selection_bin_key

    af = _tied_af()
    maf = np.where(af <= 0.5, af, 1 - af)
    he = 2 * af * (1 - af)
    n_alleles = np.full(len(af), 2)

    # `he` genuinely loses resolution: fewer distinct values than `maf` has
    assert len(np.unique(he)) < len(np.unique(maf))

    key = selection_bin_key(af, n_alleles=n_alleles, maf_k=None)
    np.testing.assert_array_equal(key, maf)


def test_a_multiallelic_site_bins_above_every_biallelic_one():
    """`maf` caps at 0.5; `1 - max(p_i)` reaches `1 - 1/k`, which is the point."""
    from plasgenomicsutils.lib.ibd_selection import selection_bin_key

    af = np.array([0.5, 0.5, 0.5])
    n_alleles = np.array([2, 3, 3])
    maf_k = np.array([0.5, 0.5, 2 / 3])          # 0.50/0.25/0.25 and 1/3 each
    key = selection_bin_key(af, n_alleles=n_alleles, maf_k=maf_k)
    assert key[0] == pytest.approx(0.5)          # biallelic: the old key, exactly
    assert key[2] == pytest.approx(2 / 3)        # more balanced than any biallelic site
    assert key[2] > key[0]


def test_the_statistic_bins_identically_on_a_biallelic_panel_with_ties():
    """End to end: the bin ids and every z-score must match, ties and all."""
    af = _tied_af(n=800, seed=3)
    he = 2 * af * (1 - af)
    mat = _matrix(n_pairs=200, n_snps=800, seed=5)
    a, abins = compute_selection_statistic(mat, af, n_bins=20)
    b, bbins = compute_selection_statistic(mat, af, n_bins=20, he=he,
                                           n_alleles=np.full(len(af), 2))
    np.testing.assert_array_equal(a["bin_id"], b["bin_id"])
    np.testing.assert_allclose(a["z_score"], b["z_score"], rtol=0, atol=0, equal_nan=True)
    assert list(abins["n_snps"]) == list(bbins["n_snps"])
