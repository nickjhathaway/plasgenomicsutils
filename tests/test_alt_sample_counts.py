"""Per-ALT carrier counts: how many *samples* carry each alternate.

"Is this variant near-private" is a question about an allele, not about a record. The
record-level reading -- count the samples that are not homozygous reference -- gives the same
answer at a biallelic site and a different one as soon as there are two alternates: a record
where every ALT is private to a different sample looks well supported, because between them
the alternates have several carriers.

`AC` cannot stand in for this. It counts alleles, so a diploidized `1/1` carrier contributes
2. And bcftools' `AC_Hom`/`AC_Het` are both zero on genuinely haploid calls, so no single
bcftools expression covers both ploidies.
"""

import pytest

pysam = pytest.importorskip("pysam")
pytest.importorskip("cyvcf2")

from plasgenomicsutils.lib.allele_counts import ALT_SAMPLE_TAG, add_alt_sample_counts

_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=c1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n"
)


def _vcf(tmp_path, body, name="in.vcf"):
    p = tmp_path / name
    p.write_text(_HDR + body)
    return str(p)


def _tag(path):
    with pysam.VariantFile(path) as vf:
        return {r.pos: tuple(r.info[ALT_SAMPLE_TAG]) for r in vf}


def test_each_alternate_is_counted_by_its_own_carriers(tmp_path):
    src = _vcf(tmp_path,
               "c1\t100\t.\tA\tC,G\t.\t.\t.\tGT\t1/1\t2/2\t0/0\t./.\n"
               "c1\t200\t.\tA\tC,G\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t./.\n"
               "c1\t300\t.\tA\tC\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t./.\n")
    out = str(tmp_path / "o.vcf")
    add_alt_sample_counts(src, out)
    assert _tag(out) == {
        100: (1, 1),    # C private to s1, G private to s2 -- two singletons
        200: (2, 0),    # C carried by two samples, G by none
        300: (2,),
    }


def test_a_het_carrier_counts_once_for_each_allele_it_names(tmp_path):
    """`1/2` is one sample carrying both alternates, not two carriers of either."""
    src = _vcf(tmp_path, "c1\t100\t.\tA\tC,G\t.\t.\t.\tGT\t1/2\t1/1\t0/0\t./.\n")
    out = str(tmp_path / "o.vcf")
    add_alt_sample_counts(src, out)
    assert _tag(out) == {100: (2, 1)}


def test_haploid_calls_are_counted_the_same_way(tmp_path):
    """The reason this is not a bcftools expression: AC_Hom and AC_Het are 0 here."""
    src = _vcf(tmp_path, "c1\t100\t.\tA\tC,G\t.\t.\t.\tGT\t1\t2\t0\t.\n")
    out = str(tmp_path / "o.vcf")
    add_alt_sample_counts(src, out)
    assert _tag(out) == {100: (1, 1)}


def test_a_biallelic_site_agrees_with_the_record_level_count(tmp_path):
    """Where the old reading was right, the new one must give the same number."""
    src = _vcf(tmp_path,
               "c1\t100\t.\tA\tC\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t./.\n"
               "c1\t200\t.\tA\tC\t.\t.\t.\tGT\t1/1\t0/0\t0/0\t./.\n")
    out = str(tmp_path / "o.vcf")
    add_alt_sample_counts(src, out)
    assert _tag(out) == {100: (2,), 200: (1,)}


def test_a_record_with_no_alt_gets_no_tag(tmp_path):
    src = _vcf(tmp_path, "c1\t100\t.\tA\t.\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t./.\n")
    out = str(tmp_path / "o.vcf")
    add_alt_sample_counts(src, out)
    with pysam.VariantFile(out) as vf:
        rec = next(iter(vf))
        assert ALT_SAMPLE_TAG not in rec.info


def test_the_singleton_filter_tests_each_alternate_separately(tmp_path):
    """The behaviour change this is all for.

    A record where C is private to s1 and G is private to s2 has two singletons and no
    well-supported allele, but the record-level count saw two non-reference samples and kept
    it. `--min-samples 1` means "drop anything carried by one sample or fewer".
    """
    from plasgenomicsutils.lib import vcf_filters as F

    src = _vcf(tmp_path,
               "c1\t100\t.\tA\tC,G\t.\t.\t.\tGT:AD\t1/1:0,9,0\t2/2:0,0,9\t"
               "0/0:9,0,0\t./.:0,0,0\n"
               "c1\t200\t.\tA\tC,G\t.\t.\t.\tGT:AD\t1/1:0,9,0\t1/1:0,9,0\t"
               "0/0:9,0,0\t./.:0,0,0\n"
               "c1\t300\t.\tA\tC\t.\t.\t.\tGT:AD\t1/1:0,9\t0/0:9,0\t"
               "0/0:9,0\t./.:0,0\n")
    out = str(tmp_path / "o.vcf")
    F.singleton_add_ads(src, out, min_samples=1)
    with pysam.VariantFile(out) as vf:
        kept = [r.pos for r in vf]
    assert kept == [200], "only the record with an alternate in two samples survives"


def test_sample_coverage_filters_singleton_recheck_is_also_per_allele(tmp_path):
    """The same test, re-run after low-coverage samples are dropped.

    It has to be recomputed rather than carried over: dropping samples changes who carries
    what, which is the whole reason this re-filter exists.
    """
    from plasgenomicsutils.lib import vcf_filters as F

    hdr = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=c1,length=100000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=ADS,Number=1,Type=Integer,Description="ADS">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n"
    )
    body = (
        # C private to s1, G private to s2 -- two singletons, no supported allele
        "c1\t100\t.\tA\tC,G\t.\t.\t.\tGT:ADS\t1/1:30\t2/2:30\t0/0:30\t0/0:30\n"
        # C carried by s1 and s2 -- supported
        "c1\t200\t.\tA\tC,G\t.\t.\t.\tGT:ADS\t1/1:30\t1/1:30\t0/0:30\t0/0:30\n"
    )
    src = str(tmp_path / "cov.vcf")
    open(src, "w").write(hdr + body)
    out = str(tmp_path / "o.vcf")
    F.sample_coverage_filter(src, out, ads_min=10, frac_min=0.5)
    with pysam.VariantFile(out) as vf:
        assert [r.pos for r in vf] == [200]


def test_the_recheck_uses_the_samples_that_remain(tmp_path):
    """Dropping the only other carrier turns a doubleton into a singleton."""
    from plasgenomicsutils.lib import vcf_filters as F

    hdr = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=c1,length=100000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=ADS,Number=1,Type=Integer,Description="ADS">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n"
    )
    # s2 carries C too, but has no coverage anywhere and will be dropped
    body = ("c1\t100\t.\tA\tC\t.\t.\t.\tGT:ADS\t1/1:30\t1/1:0\t0/0:30\t0/0:30\n"
            "c1\t200\t.\tA\tC\t.\t.\t.\tGT:ADS\t1/1:30\t1/1:0\t1/1:30\t0/0:30\n")
    src = str(tmp_path / "cov2.vcf")
    open(src, "w").write(hdr + body)
    out = str(tmp_path / "o2.vcf")
    dropped = F.sample_coverage_filter(src, out, ads_min=10, frac_min=0.5)
    assert dropped == ["s2"]
    with pysam.VariantFile(out) as vf:
        # site 100 is left with one carrier once s2 goes; site 200 still has two
        assert [r.pos for r in vf] == [200]
