"""Allele-frequency counting (global + per-region), including multiallelic + missing."""

import math

import pytest

pytest.importorskip("cyvcf2")
from plasgenomicsutils.lib.ibd_freqs import compute_allele_freqs

_VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n"
    "chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t1/1\n"      # alt counts 0,1,2
    "chr1\t20\t.\tA\tT,C\t.\t.\t.\tGT\t0/1\t1/2\t./.\n"     # multiallelic; s3 missing
)


def _af(df, snp):
    return df.set_index("snp_id")["af"].to_dict()[snp]


def test_global_and_region_allele_freqs(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    s2r = {"s1": "A", "s2": "A", "s3": "B"}
    g, r = compute_allele_freqs(str(vcf), s2r)

    # global: site10 = (0+1+2)/6 = 0.5 ; site20 = (1+2)/4 = 0.75 (s3 missing -> excluded)
    assert _af(g, "chr1:9") == pytest.approx(0.5)
    assert _af(g, "chr1:19") == pytest.approx(0.75)

    ga = r[r.group == "A"]
    gb = r[r.group == "B"]
    assert _af(ga, "chr1:9") == pytest.approx(0.25)   # s1,s2 -> (0+1)/4
    assert _af(ga, "chr1:19") == pytest.approx(0.75)   # (1+2)/4
    assert _af(gb, "chr1:9") == pytest.approx(1.0)    # s3 hom-alt
    assert math.isnan(_af(gb, "chr1:19"))              # s3 missing -> AN 0 -> NaN


def test_snp_ids_are_always_zero_based(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    g, _ = compute_allele_freqs(str(vcf))
    # VCF POS 10 and 20 -> the canonical 0-based labels; there is no toggle for 1-based
    assert list(g.snp_id) == ["chr1:9", "chr1:19"]
    assert "pos_vcf" not in g.columns              # opt-in, so the table stays small


def test_with_pos_vcf_adds_the_one_based_column(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    g, _ = compute_allele_freqs(str(vcf), with_pos_vcf=True)
    assert list(g.pos_vcf) == [10, 20]             # what the VCF itself says
    assert list(g.snp_id) == ["chr1:9", "chr1:19"]  # the key is unchanged


# ---- no metadata: frequencies over the whole file ---------------------------------
# Looking at what frequencies a callset actually holds does not need a grouping, and the
# INFO/AF field is not a substitute -- it is whatever the caller wrote, not what the
# genotypes now say after filtering and re-genotyping.

def test_without_a_grouping_only_the_whole_file_table_is_produced(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    g, r = compute_allele_freqs(str(vcf))

    assert len(g) == 2
    assert _af(g, "chr1:9") == pytest.approx(0.5)
    assert _af(g, "chr1:19") == pytest.approx(0.75)
    # the group table is empty but still shaped, so a caller can concat it either way
    assert len(r) == 0
    assert list(r.columns) == ["group", "snp_id", "n_alts", "af", "maf", "ac", "an",
                               "af_weighted", "n_samples_ad", "prevalence",
                               "n_samples_alt", "n_samples",
                               "prevalence_ad", "n_samples_alt_ad"]

    # and the whole-file numbers do not depend on whether a grouping was asked for
    with_groups, _ = compute_allele_freqs(str(vcf), {"s1": "A", "s2": "A", "s3": "B"})
    assert g.equals(with_groups)


def test_counts_travel_with_the_frequency(tmp_path):
    """af without ac/an cannot be read: 1.0 from two alleles is not 1.0 from seven hundred."""
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    g, r = compute_allele_freqs(str(vcf), {"s1": "A", "s2": "A", "s3": "B"})
    by_id = g.set_index("snp_id")

    assert by_id.loc["chr1:9", ["ac", "an"]].tolist() == [3, 6]
    # s3 is missing at site 20, so only two samples' alleles are counted
    assert by_id.loc["chr1:19", ["ac", "an"]].tolist() == [3, 4]
    assert by_id.loc["chr1:9", "af"] == pytest.approx(by_id.loc["chr1:9", "ac"]
                                                      / by_id.loc["chr1:9", "an"])
    # maf is the minor side of af
    assert by_id.loc["chr1:19", "maf"] == pytest.approx(0.25)
    assert by_id.loc["chr1:9", "maf"] == pytest.approx(0.5)

    gb = r[(r.group == "B") & (r.snp_id == "chr1:19")].iloc[0]
    assert gb.an == 0 and math.isnan(gb.af)      # nothing called: NaN, and the 0 says why


def test_the_selection_statistic_still_reads_the_widened_tables(tmp_path):
    """It selects by name, so extra columns must not disturb it."""
    from plasgenomicsutils.lib.ibd_selection import _read_af_table
    from plasgenomicsutils.lib.intervals import SNP_COORD_SYSTEM
    from plasgenomicsutils.utils.small_utils import Utils

    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)
    g, r = compute_allele_freqs(str(vcf), {"s1": "A", "s2": "A", "s3": "B"})
    gp, rp = tmp_path / "af.tsv.gz", tmp_path / "afg.tsv.gz"
    # written the way the command writes them, stamp included -- the reader rejects a
    # table with no coordinate-system stamp, which is the point of the stamp
    stamp = f"snp_coord_system={SNP_COORD_SYSTEM}"
    Utils.write_tsv_gz(g, str(gp), header_comment=stamp)
    Utils.write_tsv_gz(r, str(rp), header_comment=stamp)

    assert list(_read_af_table(str(gp), ["snp_id", "af"]).columns) == ["snp_id", "af"]
    assert list(_read_af_table(str(rp), ["group", "snp_id", "af"]).columns) == \
        ["group", "snp_id", "af"]


# ---- weighted frequency, prevalence, per-ALT --------------------------------------
# In a polyclonal infection a hard GT call throws away how much of the sample the allele
# actually is. Averaging each sample's within-sample fraction keeps it.

_MIXED = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n"
    # s1 95% alt, s2 5% alt (its 1/1 call hides that), s3 an even mix, s4 pure ref
    "chr1\t10\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:5,95\t1/1:95,5\t0/1:50,50\t0/0:100,0\n"
    # multiallelic: s1 mostly T, s2 mostly C, s3 three-way, s4 ref
    "chr1\t20\t.\tA\tT,C\t.\t.\t.\tGT:AD\t1/1:10,80,10\t2/2:10,10,80\t0/1:34,33,33\t0/0:100,0,0\n"
)


def test_weighted_af_averages_within_sample_fractions_not_hard_calls(tmp_path):
    vcf = tmp_path / "mixed.vcf"
    vcf.write_text(_MIXED)
    g, _ = compute_allele_freqs(str(vcf))
    row = g.set_index("snp_id").loc["chr1:9"]

    # hard calls: s2 is 95% REF but its 1/1 contributes two ALT alleles
    assert row.ac == 5 and row.an == 8
    assert row.af == pytest.approx(0.625)
    # weighted: (0.95 + 0.05 + 0.50 + 0.00) / 4, over samples rather than alleles
    assert row.af_weighted == pytest.approx(0.375)
    assert row.n_samples_ad == 4
    # a deep sample does not outvote a shallow one: each contributes at most 1
    assert 0.0 <= row.af_weighted <= 1.0


def test_prevalence_counts_samples_carrying_the_allele(tmp_path):
    vcf = tmp_path / "mixed.vcf"
    vcf.write_text(_MIXED)
    g, _ = compute_allele_freqs(str(vcf))
    row = g.set_index("snp_id").loc["chr1:9"]
    assert row.n_samples_alt == 3 and row.n_samples == 4   # s1, s2, s3 carry it
    assert row.prevalence == pytest.approx(0.75)
    # sample counts are not allele counts -- at ploidy 2 they differ by a factor of two
    assert row.an == 2 * row.n_samples


def test_per_alt_splits_a_multiallelic_site_consistently(tmp_path):
    vcf = tmp_path / "mixed.vcf"
    vcf.write_text(_MIXED)
    coll, _ = compute_allele_freqs(str(vcf))
    per, _ = compute_allele_freqs(str(vcf), per_alt=True)

    multi = per[per.snp_id == "chr1:19"]
    assert multi.alt.tolist() == ["T", "C"]
    assert multi.alt_index.tolist() == [1, 2]
    assert multi.an.nunique() == 1                     # the denominator is the site's

    c = coll.set_index("snp_id").loc["chr1:19"]
    assert multi.ac.sum() == c.ac                      # the ALTs partition the alt count
    assert multi.af_weighted.sum() == pytest.approx(c.af_weighted)
    # (0.80 + 0.10 + 0.33 + 0.00) / 4 and (0.10 + 0.80 + 0.33 + 0.00) / 4
    assert multi.af_weighted.tolist() == pytest.approx([0.3075, 0.3075])
    assert multi.n_samples_alt.tolist() == [2, 1]      # GT-based: T in s1+s3, C in s2

    # a biallelic site reads the same either way
    assert per[per.snp_id == "chr1:9"].af.iloc[0] == pytest.approx(
        coll.set_index("snp_id").loc["chr1:9"].af)


def test_weighted_is_nan_without_ad_and_can_be_turned_off(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)                                # no FORMAT/AD at all
    g, _ = compute_allele_freqs(str(vcf))
    assert g.af_weighted.isna().all()
    assert (g.n_samples_ad == 0).all()
    assert g.af.notna().any()                           # the count-based columns still work

    vcf2 = tmp_path / "mixed.vcf"
    vcf2.write_text(_MIXED)
    off, _ = compute_allele_freqs(str(vcf2), weighted=False)
    assert off.af_weighted.isna().all()
    on, _ = compute_allele_freqs(str(vcf2))
    assert on.af_weighted.notna().all()
    # prevalence_ad reads the same AD, so it goes too -- everything else is untouched
    assert off.prevalence_ad.isna().all()
    ad_derived = ("af_weighted", "n_samples_ad", "prevalence_ad", "n_samples_alt_ad")
    same = [c for c in on.columns if c not in ad_derived]
    assert on[same].equals(off[same])


def test_group_tables_carry_the_same_columns(tmp_path):
    vcf = tmp_path / "mixed.vcf"
    vcf.write_text(_MIXED)
    g, r = compute_allele_freqs(str(vcf), {"s1": "A", "s2": "A", "s3": "B", "s4": "B"},
                                per_alt=True)
    assert set(g.columns).issubset(set(r.columns))
    a = r[(r.group == "A") & (r.snp_id == "chr1:9")].iloc[0]
    # group A is s1 (95% alt) and s2 (5% alt): weighted 0.5, but both called 1/1 -> af 1.0
    assert a.af == pytest.approx(1.0)
    assert a.af_weighted == pytest.approx(0.5)
    assert a.prevalence == pytest.approx(1.0) and a.n_samples == 2


# ---- AD-based prevalence ----------------------------------------------------------
# A minor clone can be well supported by reads and still absent from the genotype, because
# the caller committed to one. `prevalence` and `prevalence_ad` are the two readings.

_CLONE = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\ts5\n"
    # s1 hom-alt; s2/s3 called 0/0 but carry a 3% and 5% clone; s4 has one stray alt
    # read (below the 2-read floor); s5 is pure reference
    "chr1\t30\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,100\t0/0:97,3\t0/0:95,5\t0/0:99,1\t0/0:100,0\n"
)


def test_ad_prevalence_finds_clones_the_genotype_missed(tmp_path):
    vcf = tmp_path / "clone.vcf"
    vcf.write_text(_CLONE)
    r = compute_allele_freqs(str(vcf))[0].iloc[0]

    assert r.n_samples_alt == 1 and r.prevalence == pytest.approx(0.2)      # GT: s1 only
    assert r.n_samples_alt_ad == 3 and r.prevalence_ad == pytest.approx(0.6)  # + s2, s3
    # s4's single alt read is below the 2-read floor, so it is not "carrying" it
    assert r.n_samples_ad == 5


def test_the_ad_presence_thresholds_are_tunable(tmp_path):
    vcf = tmp_path / "clone.vcf"
    vcf.write_text(_CLONE)

    def prev(**kw):
        return compute_allele_freqs(str(vcf), **kw)[0].iloc[0].prevalence_ad

    assert prev() == pytest.approx(0.6)                        # 3% and 5% both count
    assert prev(ad_min_freq=0.04) == pytest.approx(0.4)        # 3% drops out
    assert prev(ad_min_freq=0.10) == pytest.approx(0.2)        # both drop out
    # the read floor is the other half: at 6 reads only s1 (100 reads) survives
    assert prev(ad_min_reads=6) == pytest.approx(0.2)
    # both floors have to be met, not either
    assert prev(ad_min_reads=4, ad_min_freq=0.01) == pytest.approx(0.4)


def test_ad_prevalence_is_nan_without_ad(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(_VCF)                                       # no FORMAT/AD
    g, _ = compute_allele_freqs(str(vcf))
    assert g.prevalence_ad.isna().all()
    assert (g.n_samples_alt_ad == 0).all()
    assert g.prevalence.notna().any()                          # the GT reading still works


def test_ad_prevalence_is_per_alt_when_asked(tmp_path):
    vcf = tmp_path / "mixed.vcf"
    vcf.write_text(_MIXED)
    per, _ = compute_allele_freqs(str(vcf), per_alt=True)
    multi = per[per.snp_id == "chr1:19"]
    # ALT=T is >=1% of s1 (80%), s2 (10%) and s3 (33%); ALT=C likewise
    assert multi.n_samples_alt_ad.tolist() == [3, 3]
    # but the genotypes name only one or two of them
    assert multi.n_samples_alt.tolist() == [2, 1]


# ---- n_alts ----------------------------------------------------------------------

_MULTI3 = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n"
    "chr1\t10\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:5,95\t1/1:95,5\t0/1:50,50\t0/0:100,0\n"
    "chr1\t20\t.\tA\tT,C\t.\t.\t.\tGT:AD\t1/1:10,80,10\t2/2:10,10,80\t0/1:34,33,33\t0/0:100,0,0\n"
    "chr1\t30\t.\tA\tT,C,G\t.\t.\t.\tGT:AD\t1/1:10,70,10,10\t2/2:10,10,70,10\t3/3:10,10,10,70\t0/0:100,0,0,0\n"
)


def test_n_alts_flags_a_multiallelic_site_in_the_collapsed_table(tmp_path):
    vcf = tmp_path / "multi3.vcf"
    vcf.write_text(_MULTI3)
    g, _ = compute_allele_freqs(vcf.as_posix())

    assert g.set_index("snp_id")["n_alts"].to_dict() == {
        "chr1:9": 1, "chr1:19": 2, "chr1:29": 3}
    # without it these sites are hard to tell apart: every one is a 0.75 prevalence
    assert g.prevalence.nunique() == 1


def test_n_alts_is_the_site_total_on_every_per_alt_row(tmp_path):
    vcf = tmp_path / "multi3.vcf"
    vcf.write_text(_MULTI3)
    per, _ = compute_allele_freqs(vcf.as_posix(), per_alt=True)

    # it names how many rows the site has, so "1 of 3" is readable from one row
    assert (per.groupby("snp_id").size() == per.groupby("snp_id").n_alts.first()).all()
    assert per[per.snp_id == "chr1:29"].alt.tolist() == ["T", "C", "G"]
    assert per[per.snp_id == "chr1:29"].alt_index.tolist() == [1, 2, 3]

    # the split still partitions the collapsed counts at three ALTs
    coll, _ = compute_allele_freqs(vcf.as_posix())
    assert (per.groupby("snp_id").ac.sum().sort_index()
            == coll.set_index("snp_id").ac.sort_index()).all()


def test_n_alts_is_carried_by_the_group_table_too(tmp_path):
    vcf = tmp_path / "multi3.vcf"
    vcf.write_text(_MULTI3)
    _, r = compute_allele_freqs(vcf.as_posix(),
                                {"s1": "A", "s2": "A", "s3": "B", "s4": "B"})
    assert "n_alts" in r.columns
    # a site-level property, so it does not vary by group
    assert (r.groupby("snp_id").n_alts.nunique() == 1).all()
