"""wsaf_profile: does a sample have a dominant clone, and what filter reduces it to one?"""

import numpy as np
import pytest

pytest.importorskip("cyvcf2")
from plasgenomicsutils.lib.wsaf import (  # noqa: E402
    WSAF_BIN,
    _count_bands,
    _mode_of,
    _quantile_of,
    classify_profile,
    min_freq_needed,
    rate_at_or_above,
    wsaf_profile,
)

HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=4000000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
]


def _write_vcf(path, samples, rows):
    """rows: list of (pos, [(ref_ad, alt_ad), ...] per sample)."""
    lines = list(HEADER)
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
    for pos, ads in rows:
        cells = []
        for ref, alt in ads:
            gt = "0/1" if ref and alt else ("1/1" if alt else "0/0")
            cells.append(f"{gt}:{ref},{alt}")
        lines.append(f"chr1\t{pos}\t.\tA\tT\t100\t.\t.\tGT:AD\t" + "\t".join(cells))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _mixture(rng, minor_prop, n_sites, depth=200, het_share=0.15):
    """A sample whose het sites carry the minor strain at `minor_prop` of the reads."""
    out = []
    for _ in range(n_sites):
        if rng.random() < het_share:
            alt = rng.binomial(depth, minor_prop)
            out.append((depth - alt, alt))
        else:
            out.append((depth, 0))
    return out


def _hist(fracs):
    """The per-sample histogram every summary is derived from."""
    n = int(round(0.5 / WSAF_BIN))
    return np.histogram(fracs, bins=np.arange(0, 0.5 + WSAF_BIN / 2, WSAF_BIN))[0][:n]


# ---- the binned summaries ---------------------------------------------------------

def test_mode_and_quantiles_come_off_the_histogram():
    rng = np.random.default_rng(0)
    v = np.concatenate([rng.normal(0.20, 0.01, 500), rng.uniform(0, 0.05, 50)])
    h = _hist(v)
    assert _mode_of(h) == pytest.approx(0.20, abs=WSAF_BIN)
    assert _quantile_of(h, 0.5) == pytest.approx(np.median(v), abs=2 * WSAF_BIN)
    assert np.isnan(_mode_of(np.zeros(50)))
    assert np.isnan(_quantile_of(np.zeros(50), 0.5))


def test_count_bands_separates_peaks():
    rng = np.random.default_rng(1)
    one = rng.normal(0.10, 0.005, 400)
    two = np.concatenate([one, rng.normal(0.40, 0.005, 400)])
    assert _count_bands(_hist(one)) == 1
    assert _count_bands(_hist(two)) == 2
    assert _count_bands(np.zeros(50)) == 0


def test_rate_at_or_above_is_over_covered_sites():
    # 30 het sites at 0.40 in a 1000-site panel
    h = _hist(np.full(30, 0.40))
    assert rate_at_or_above(h, 1000, 0.35) == pytest.approx(0.03)
    assert rate_at_or_above(h, 1000, 0.45) == pytest.approx(0.0)
    assert np.isnan(rate_at_or_above(h, 0, 0.35))


# ---- the filter the sample needs --------------------------------------------------

def test_min_freq_needed_is_the_threshold_that_clears_the_residue():
    # 100 sites at 0.15 in a 1000-site panel: a filter at 0.15 leaves 10% of sites het,
    # one just above it leaves none
    h = _hist(np.full(100, 0.15))
    assert min_freq_needed(h, 1000, max_residual_het=0.02) == pytest.approx(0.16, abs=1e-9)
    # tolerate that much residue and no filtering is needed at all
    assert min_freq_needed(h, 1000, max_residual_het=0.20) == pytest.approx(0.02, abs=1e-9)


def test_min_freq_needed_is_nan_for_a_co_dominant_strain():
    # a band at 0.49 cannot be filtered away: a threshold above it deletes the dominant call
    h = _hist(np.full(500, 0.49))
    assert np.isnan(min_freq_needed(h, 1000, max_residual_het=0.02))


def test_a_clonal_sample_needs_no_filter():
    # a scattering of noise-level het calls, well under the residue bar
    h = _hist(np.full(5, 0.30))
    assert min_freq_needed(h, 20000, max_residual_het=0.02) == pytest.approx(0.02, abs=1e-9)


# ---- the classifier --------------------------------------------------------------

def test_classify_names_what_it_would_take():
    assert classify_profile(0.02, 20000) == "monoclonal"
    assert classify_profile(0.20, 20000) == "dominant_clone"      # <= 1 - 0.70
    assert classify_profile(0.30, 20000) == "dominant_clone"      # exactly at the bar
    assert classify_profile(0.35, 20000) == "mixed"               # a harder filter than asked
    assert classify_profile(float("nan"), 20000) == "mixed"       # no filter works at all
    assert classify_profile(0.02, 200) == "undetermined"          # thin coverage


def test_min_dominant_is_the_callers_choice():
    # the same sample, judged against two compositions
    assert classify_profile(0.35, 20000, min_dominant=0.70) == "mixed"
    assert classify_profile(0.35, 20000, min_dominant=0.60) == "dominant_clone"


# ---- end to end over a file ------------------------------------------------------

def test_dominant_frac_recovers_the_strain_proportions(tmp_path):
    rng = np.random.default_rng(7)
    props = {"dom_10": 0.10, "dom_25": 0.25, "even": 0.49}
    per = {s: _mixture(rng, p, 800) for s, p in props.items()}
    samples = list(props)
    rows = [(1000 + 100 * i, [per[s][i] for s in samples]) for i in range(800)]
    vcf = _write_vcf(tmp_path / "mix.vcf", samples, rows)

    df = wsaf_profile(vcf, min_sites=100).set_index("sample")
    for s, p in props.items():
        assert df.loc[s, "minor_mode"] == pytest.approx(p, abs=0.03), s
    # the filter each one needs sits just above its minor strain's proportion, because the
    # binomial spread around it has to be cleared too
    for s, p in (("dom_10", 0.10), ("dom_25", 0.25)):
        assert p < df.loc[s, "min_freq_needed"] <= p + 0.08, s
        assert df.loc[s, "class"] == "dominant_clone", s
    # an even mixture cannot be filtered at all
    assert np.isnan(df.loc["even", "min_freq_needed"])
    assert df.loc["even", "class"] == "mixed"
    assert df.loc["dom_10", "dominant_frac"] > df.loc["dom_25", "dominant_frac"]


def test_residual_het_rate_predicts_what_the_filter_leaves(tmp_path):
    # a quarter of sites het at 0.40, which a 0.30 filter does NOT remove
    rows = [(1000 + i, [(120, 80) if i % 4 == 0 else (200, 0)]) for i in range(400)]
    vcf = _write_vcf(tmp_path / "resid.vcf", ["s"], rows)
    df = wsaf_profile(vcf, min_sites=100, min_dominant=0.70)
    assert df["residual_het_rate"].iloc[0] == pytest.approx(0.25, abs=0.01)
    assert df["class"].iloc[0] == "mixed"


def test_minor_fraction_stays_within_half_at_a_multiallelic_site(tmp_path):
    # three alleles at a third each: the runner-up's share is 0.33, not 0.67
    lines = list(HEADER)
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1")
    lines.append("chr1\t1000\t.\tA\tT,G\t100\t.\t.\tGT:AD\t0/1:30,30,30")
    p = tmp_path / "multi.vcf"
    p.write_text("\n".join(lines) + "\n")

    sites = tmp_path / "sites.tsv.gz"
    wsaf_profile(str(p), min_sites=1, sites_out=str(sites))
    import pandas as pd
    got = pd.read_csv(sites, sep="\t")
    assert len(got) == 1
    assert got["minor_frac"].iloc[0] == pytest.approx(1 / 3, abs=1e-4)


def test_low_depth_sites_are_not_counted(tmp_path):
    # at 6x the fractions are multiples of 1/6, too coarse to place a band
    rows = [(1000, [(3, 3)]), (2000, [(50, 50)])]
    vcf = _write_vcf(tmp_path / "depth.vcf", ["s1"], rows)
    df = wsaf_profile(vcf, min_depth=10, min_sites=1)
    assert df["n_sites"].iloc[0] == 1
    assert df["n_het"].iloc[0] == 1


def test_a_single_stray_read_is_not_a_heterozygous_site(tmp_path):
    # 1 read of 10 is 10% -- over min_minor, but one read is what error produces
    rows = [(1000 + i, [(9, 1)]) for i in range(200)]
    vcf = _write_vcf(tmp_path / "stray.vcf", ["s1"], rows)
    assert wsaf_profile(vcf, min_sites=1)["n_het"].iloc[0] == 0
    assert wsaf_profile(vcf, min_sites=1, min_minor_reads=1)["n_het"].iloc[0] == 200


def test_wsmaf_flips_with_the_population_frequency(tmp_path):
    # the ALT is the population-level MAJOR allele here (3 of 4 samples are hom-alt),
    # so at the het sample WSMAF is the REF fraction, not the ALT fraction
    samples = ["het", "a", "b", "c"]
    rows = [(1000 + i, [(80, 20), (0, 100), (0, 100), (0, 100)]) for i in range(200)]
    vcf = _write_vcf(tmp_path / "flip.vcf", samples, rows)
    sites = tmp_path / "sites.tsv.gz"
    df = wsaf_profile(vcf, min_sites=1, sites_out=str(sites)).set_index("sample")

    import pandas as pd
    got = pd.read_csv(sites, sep="\t")
    got = got[got["sample"] == "het"]
    assert got["alt_frac"].mean() == pytest.approx(0.20, abs=0.01)
    assert (got["plaf"] > 0.5).all()
    assert got["wsmaf"].mean() == pytest.approx(0.80, abs=0.01)
    assert df.loc["het", "wsmaf_mean"] == pytest.approx(0.80, abs=0.01)


def test_sites_out_agrees_with_the_summary(tmp_path):
    rng = np.random.default_rng(11)
    samples = ["a", "b"]
    per = {s: _mixture(rng, 0.15, 400) for s in samples}
    rows = [(1000 + i, [per[s][i] for s in samples]) for i in range(400)]
    vcf = _write_vcf(tmp_path / "pair.vcf", samples, rows)
    sites = tmp_path / "sites.tsv.gz"
    df = wsaf_profile(vcf, min_sites=10, sites_out=str(sites)).set_index("sample")

    import pandas as pd
    long = pd.read_csv(sites, sep="\t")
    assert set(long.columns) == {"sample", "snp_id", "minor_frac", "alt_frac", "plaf", "wsmaf"}
    for s in samples:
        m = long.loc[long["sample"] == s, "minor_frac"]
        assert len(m) == df.loc[s, "n_het"]
        assert m.median() == pytest.approx(df.loc[s, "minor_median"], abs=2 * WSAF_BIN)
    # snp_id is zero-based, the convention used everywhere else in the package, so the
    # first record's VCF POS of 1000 is written as 999
    assert long["snp_id"].iloc[0] == "chr1:999"


def test_a_clonal_sample_is_called_monoclonal(tmp_path):
    # every site hom-ref except a handful of stray reads
    rows = [(1000 + i, [(100, 4) if i % 400 == 0 else (100, 0)]) for i in range(2000)]
    vcf = _write_vcf(tmp_path / "clonal.vcf", ["s1"], rows)
    df = wsaf_profile(vcf)
    assert df["n_het"].iloc[0] > 0            # not zero: it is the rate that is low
    assert df["class"].iloc[0] == "monoclonal"
    assert df["min_freq_needed"].iloc[0] == pytest.approx(0.02, abs=1e-9)
