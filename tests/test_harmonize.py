"""Tests for harmonize duplicate-position handling and ALT-union logic."""

import pytest

pysam = pytest.importorskip("pysam")

from plasgenomicsutils.lib import harmonize as H


VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)


def _write_vcf(path, body):
    path.write_text(VCF_HEADER + body)
    return str(path)


def _records(path):
    with pysam.VariantFile(path) as vf:
        return [(r.chrom, r.pos, r.alleles, {s: tuple(r.samples[s]["GT"]) for s in r.samples})
                for r in vf]


def test_duplicate_snp_and_noalt_keeps_the_snp(tmp_path):
    # A real SNP record and an overlapping no-ALT record at the same position.
    f = _write_vcf(tmp_path / "a.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/1:6,4\n"
        "chr1\t100\t.\tA\t.\t.\t.\t.\tGT:AD\t0/0:9\t0/0:8\n"
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,9\t0/0:10,0\n")

    union, dups, amb = H.accumulate_union([f], min_ad=3, min_af=0.005, het_min_af=0.2)
    assert union[("chr1", 100)] == ["A", "T"]     # SNP's ALT survives
    assert len(dups) == 1                          # the duplicate was noticed
    assert len(amb) == 0                           # not ambiguous (only one real-ALT record)

    out = str(tmp_path / "out.vcf")
    H.harmonize_file(f, out, union, 3, 0.005, 0.2)
    recs = _records(out)
    at_100 = [r for r in recs if r[1] == 100]
    assert len(at_100) == 1                        # exactly one record at the position
    _, _, alleles, gts = at_100[0]
    assert alleles == ("A", "T")
    assert gts["S1"][:2] == (0, 1)                 # the real SNP genotype is kept, not 0/0


def test_two_real_alt_records_flagged_ambiguous(tmp_path):
    f = _write_vcf(tmp_path / "b.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"
        "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD\t1/1:0,8\t0/1:4,4\n")
    union, dups, amb = H.accumulate_union([f], min_ad=3, min_af=0.005, het_min_af=0.2)
    assert len(dups) == 1
    assert len(amb) == 1                           # two ALT-bearing records -> ambiguous


def test_multiallelic_reduction_keeps_ad_length_consistent(tmp_path):
    # G>T,A where A has zero support (dropped) and one sample has missing AD.
    # After reduction to G>T every sample's AD must have exactly 2 values, or
    # bcftools merge fails with "Incorrect number of FORMAT/AD values".
    f = _write_vcf(tmp_path / "m.vcf",
        "chr1\t100\t.\tG\tT,A\t.\t.\t.\tGT:AD\t0/1:10,5,0\t./.:.\t1/1:0,9,0\n")
    union, _dups, _amb = H.accumulate_union([f], min_ad=3, min_af=0.005, het_min_af=0.2)
    assert union[("chr1", 100)] == ["G", "T"]        # A dropped (no support)

    out = str(tmp_path / "out.vcf")
    H.harmonize_file(f, out, union, 3, 0.005, 0.2)
    with pysam.VariantFile(out) as vf:
        rec = next(iter(vf))
        assert len(rec.alleles) == 2
        for s in rec.samples:
            ad = rec.samples[s]["AD"]
            assert len(ad) == 2                      # consistent Number=R length


def test_indel_context_records_dropped_by_default(tmp_path):
    # A SNP plus a no-ALT INDEL-flagged record and a multi-base REF indel.
    body = (
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"
        "chr1\t100\t.\tA\t.\t.\t.\tINDEL\tGT:AD\t0/0:9\t0/0:8\n"      # no-ALT indel: bcftools misses it
        "chr1\t200\t.\tATG\t.\t.\t.\tINDEL\tGT:AD\t0/0:7\t0/0:6\n")   # multi-base REF indel
    f = _write_vcf(tmp_path / "i.vcf", body)

    # default: indels dropped -> no duplicate at 100, no site at 200
    union, dups, _amb = H.accumulate_union([f], 3, 0.005, 0.2)
    assert union[("chr1", 100)] == ["A", "T"]
    assert ("chr1", 200) not in union
    assert len(dups) == 0                            # the indel record never entered

    # is_indel_context flags them
    with pysam.VariantFile(f) as vf:
        recs = list(vf)
    assert H.is_indel_context(recs[1]) is True       # no-ALT INDEL flag
    assert H.is_indel_context(recs[2]) is True       # multi-base REF
    assert H.is_indel_context(recs[0]) is False      # the SNP


def test_stale_format_fields_detection(tmp_path):
    # PL (Number=G) is stale after allele reshaping; AD (R) is maintained; GT/DP are not per-allele.
    hdr = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD:DP:PL\t0/1:5,5:10:120,0,130\n"
    )
    f = tmp_path / "pl.vcf"
    f.write_text(hdr)
    assert H.stale_format_fields(str(f)) == ["PL"]


def test_alt_union_across_files(tmp_path):
    a = _write_vcf(tmp_path / "a.vcf", "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t1/1:0,9\n")
    b = _write_vcf(tmp_path / "b.vcf", "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD\t1/1:0,7\t0/0:8,0\n")
    union, _dups, _amb = H.accumulate_union([a, b], min_ad=3, min_af=0.005, het_min_af=0.2)
    assert union[("chr1", 100)] == ["A", "G", "T"]  # union of ALTs across files, sorted
