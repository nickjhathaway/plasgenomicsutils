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

    union, dups, amb, _st = H.accumulate_union([f], min_ad=3, min_af=0.005, het_min_af=0.2)
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
    union, dups, amb, _st = H.accumulate_union([f], min_ad=3, min_af=0.005, het_min_af=0.2)
    assert len(dups) == 1
    assert len(amb) == 1                           # two ALT-bearing records -> ambiguous


def test_multiallelic_reduction_keeps_ad_length_consistent(tmp_path):
    # G>T,A where A has zero support (dropped) and one sample has missing AD.
    # After reduction to G>T every sample's AD must have exactly 2 values, or
    # bcftools merge fails with "Incorrect number of FORMAT/AD values".
    f = _write_vcf(tmp_path / "m.vcf",
        "chr1\t100\t.\tG\tT,A\t.\t.\t.\tGT:AD\t0/1:10,5,0\t./.:.\t1/1:0,9,0\n")
    union, _dups, _amb, _st = H.accumulate_union([f], min_ad=3, min_af=0.005, het_min_af=0.2)
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
    union, dups, _amb, _st = H.accumulate_union([f], 3, 0.005, 0.2)
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
    union, _dups, _amb, _st = H.accumulate_union([a, b], min_ad=3, min_af=0.005, het_min_af=0.2)
    assert union[("chr1", 100)] == ["A", "G", "T"]  # union of ALTs across files, sorted


# --------------------------------------------------------------------------- #
#  Reporting counters                                                          #
# --------------------------------------------------------------------------- #


def test_pass1_counts_what_cleaning_removed(tmp_path):
    """The cleaning tally is how you tell a threshold doing useful work from one quietly
    discarding real alleles, so the counts have to be real, not derived from the output."""
    # site 100: a 1-read ALT in one sample only -> zeroed, record becomes ref-only
    # site 200: a solid ALT -> survives
    # site 300: two ALTs, one solid one spurious -> one removed, record keeps an ALT
    f = _write_vcf(tmp_path / "a.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:40,1\t0/0:38,0\n"
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,30\t0/0:25,0\n"
        "chr1\t300\t.\tA\tT,G\t.\t.\t.\tGT:AD\t1/1:0,30,1\t0/0:25,0,0\n")

    union, _d, _a, st = H.accumulate_union([f], min_ad=3, min_af=0.005, het_min_af=0.2)
    s = st["per_file"][f]
    assert s["processed"] == 3
    assert s["sites"] == 3
    assert s["alts_removed"] == 2                 # the 1-read ALT at 100 and the one at 300
    assert s["reduced_to_ref_only"] == 1          # site 100 only
    assert st["union_sites"] == 3
    assert st["union_with_alts"] == 2             # 200 and 300
    assert st["union_dropped"] == 1               # 100 has no real ALT anywhere


def test_pass2_counts_written_added_and_absent(tmp_path):
    """`absent` is the count that explains a later merge full of missing AD, so it has to be
    the union sites this file holds no record for -- not a difference of two other numbers."""
    a = _write_vcf(tmp_path / "a.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,30\t0/0:25,0\n"
        "chr1\t200\t.\tA\tG\t.\t.\t.\tGT:AD\t1/1:0,28\t0/0:22,0\n")
    b = _write_vcf(tmp_path / "b.vcf",
        "chr1\t100\t.\tA\tC\t.\t.\t.\tGT:AD\t1/1:0,26\t0/0:24,0\n"
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,27\t0/0:21,0\n")

    union, _d, _a, _st = H.accumulate_union([a, b], min_ad=3, min_af=0.005, het_min_af=0.2)
    assert set(union) == {("chr1", 100), ("chr1", 200), ("chr1", 300)}
    assert union[("chr1", 100)] == ["A", "C", "T"]        # both ALTs unioned

    sa = H.harmonize_file(a, str(tmp_path / "a.out.vcf"), union, 3, 0.005, 0.2)
    assert sa["written"] == 2
    assert sa["alts_added"] == 1        # site 100 gained C
    assert sa["absent"] == 1           # 300 is not in this file
    sb = H.harmonize_file(b, str(tmp_path / "b.out.vcf"), union, 3, 0.005, 0.2)
    assert sb["written"] == 2
    assert sb["alts_added"] == 1        # site 100 gained T
    assert sb["absent"] == 1           # 200 is not in this file


def test_absent_sites_become_missing_AD_after_a_merge(tmp_path):
    """Why `absent` is worth printing: the merge fills those samples with missing FORMAT/AD,
    which is what makes an integer-AD reader (hmmibd-rs) fail."""
    import shutil
    import subprocess

    if shutil.which("bcftools") is None:
        pytest.skip("bcftools not on PATH")

    a = _write_vcf(tmp_path / "a.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,30\t0/0:25,0\n")
    b = tmp_path / "b.vcf"
    b.write_text(VCF_HEADER.replace("S1\tS2", "S3\tS4") +
                 "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,26\t0/0:24,0\n"
                 "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,27\t0/0:21,0\n")

    union, _d, _a, _st = H.accumulate_union([a, str(b)], 3, 0.005, 0.2)
    sa = H.harmonize_file(a, str(tmp_path / "ha.vcf"), union, 3, 0.005, 0.2)
    H.harmonize_file(str(b), str(tmp_path / "hb.vcf"), union, 3, 0.005, 0.2)
    assert sa["absent"] == 1                      # file a has no record at 300

    for n in ("ha", "hb"):
        subprocess.run(f"bcftools view {tmp_path}/{n}.vcf -Oz -o {tmp_path}/{n}.vcf.gz "
                       f"&& bcftools index -f {tmp_path}/{n}.vcf.gz",
                       shell=True, check=True, executable="/bin/bash")
    merged = subprocess.run(
        f"bcftools merge {tmp_path}/ha.vcf.gz {tmp_path}/hb.vcf.gz | "
        f"bcftools query -f '%POS[\\t%AD]\\n'",
        shell=True, capture_output=True, text=True, executable="/bin/bash").stdout

    at_300 = [l for l in merged.splitlines() if l.startswith("300")][0]
    assert "\t." in at_300                        # file a's samples have missing AD there
    # ... and the documented filter removes exactly that record
    kept = subprocess.run(
        f"bcftools merge {tmp_path}/ha.vcf.gz {tmp_path}/hb.vcf.gz | "
        f"bcftools view -H -e 'FMT/AD=\".\"' | wc -l",
        shell=True, capture_output=True, text=True, executable="/bin/bash").stdout
    assert int(kept.strip()) == 1                 # only site 100 survives
