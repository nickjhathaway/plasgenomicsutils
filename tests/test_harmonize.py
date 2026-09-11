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


# ---- the per-allele fields move with the alleles -----------------------------------

R_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##INFO=<ID=AD,Number=R,Type=Integer,Description="x">\n'
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="x">\n'
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="x">\n'
    '##INFO=<ID=RPBZ,Number=1,Type=Float,Description="x">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    '##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="ADF">\n'
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)


def _write_r(path, body, samples="S1\tS2"):
    path.write_text(R_HEADER.replace("S1\tS2", samples) + body)
    return str(path)


def _first(path):
    with pysam.VariantFile(path) as vf:
        return [(r.alleles, dict(r.info), {s: (tuple(r.samples[s]["GT"]), r.samples[s].phased,
                                             tuple(r.samples[s]["AD"]),
                                             tuple(r.samples[s]["ADF"]))
                                        for s in r.samples}) for r in vf]


def test_every_number_r_field_is_relaid_on_the_union_and_gt_is_reindexed(tmp_path):
    """Two files that saw different ALTs at one site. The union is A>C,T; file a's T moves
    to slot 2 in FORMAT/AD, FORMAT/ADF *and* INFO/AD, C gets a zero, and the genotypes
    are re-indexed (1/1 for T -> 2/2) rather than re-called, phasing intact. INFO/AC
    (Number=A) cannot be re-laid and goes."""
    a = _write_r(tmp_path / "a.vcf",
                 "chr1\t100\t.\tA\tT\t.\t.\tAD=25,30;AC=2;DP=55;RPBZ=1.5\tGT:AD:ADF:PL\t"
                 "1|1:0,30:0,16:90,9,0\t0/0:25,0:12,0:0,9,90\n")
    b = _write_r(tmp_path / "b.vcf",
                 "chr1\t100\t.\tA\tC\t.\t.\tAD=24,26;AC=2;DP=50;RPBZ=-0.5\tGT:AD:ADF:PL\t"
                 "1/1:0,26:0,13:90,9,0\t0/0:24,0:11,0:0,9,90\n", samples="S3\tS4")
    union, _d, _a, _st = H.accumulate_union([a, b], 0, 0.0, 0.2)
    assert union[("chr1", 100)] == ["A", "C", "T"]

    H.harmonize_file(a, str(tmp_path / "ha.vcf"), union, 0, 0.0, 0.2,
                     regenotype=False, stale_info=("AC",))
    (alleles, info, samples), = _first(str(tmp_path / "ha.vcf"))
    assert alleles == ("A", "C", "T")
    assert tuple(info["AD"]) == (25, 0, 30) and "AC" not in info
    assert info["DP"] == 55 and abs(info["RPBZ"] - 1.5) < 1e-6    # Number=1 stays
    assert samples["S1"] == ((2, 2), True, (0, 0, 30), (0, 0, 16))
    assert samples["S2"] == ((0, 0), False, (25, 0, 0), (12, 0, 0))

    H.harmonize_file(b, str(tmp_path / "hb.vcf"), union, 0, 0.0, 0.2,
                     regenotype=False, stale_info=("AC",))
    (alleles, info, samples), = _first(str(tmp_path / "hb.vcf"))
    assert tuple(info["AD"]) == (24, 26, 0)
    assert samples["S3"] == ((1, 1), False, (0, 26, 0), (0, 13, 0))


def test_regenotype_false_keeps_the_callers_genotype(tmp_path):
    """A 27,3 het: the caller said 0/1 from likelihoods. Re-genotyping from AD at
    het_min_af=0.2 would make it 0/0; regenotype=False re-indexes it and leaves it be."""
    a = _write_r(tmp_path / "a.vcf",
                 "chr1\t100\t.\tA\tT\t.\t.\tAD=27,3;DP=30\tGT:AD:ADF\t0/1:27,3:14,2\t0/0:1,0:1,0\n")
    b = _write_r(tmp_path / "b.vcf",
                 "chr1\t100\t.\tA\tC\t.\t.\tAD=0,10;DP=10\tGT:AD:ADF\t1/1:0,10:0,5\t0/0:1,0:1,0\n",
                 samples="S3\tS4")
    union, _d, _a, _st = H.accumulate_union([a, b], 0, 0.0, 0.2)

    H.harmonize_file(a, str(tmp_path / "recalled.vcf"), union, 0, 0.0, 0.2)
    (_al, _info, recalled), = _first(str(tmp_path / "recalled.vcf"))
    assert recalled["S1"][0] == (0, 0)                        # re-called from AD

    H.harmonize_file(a, str(tmp_path / "kept.vcf"), union, 0, 0.0, 0.2, regenotype=False)
    (_al, _info, kept), = _first(str(tmp_path / "kept.vcf"))
    assert kept["S1"][0] == (0, 2)                            # 0/1 for T, now slot 2


def test_the_same_alts_in_another_order_still_reindex_the_genotypes(tmp_path):
    """A file listing T,C against a union of C,T gained nothing, but its 0/1 meant T and
    must come out 0/2 -- whichever way genotypes are otherwise handled."""
    a = _write_r(tmp_path / "a.vcf",
                 "chr1\t100\t.\tA\tT,C\t.\t.\tAD=20,10,5;DP=35\tGT:AD:ADF\t"
                 "0/1:10,10,0:5,5,0\t0/2:10,0,5:5,0,3\n")
    b = _write_r(tmp_path / "b.vcf",
                 "chr1\t100\t.\tA\tC,T\t.\t.\tAD=20,5,10;DP=35\tGT:AD:ADF\t"
                 "0/1:10,5,0:5,3,0\t0/2:10,0,10:5,0,5\n", samples="S3\tS4")
    union, _d, _a, _st = H.accumulate_union([a, b], 0, 0.0, 0.2)
    assert union[("chr1", 100)] == ["A", "C", "T"]
    for regenotype in (True, False):
        out = str(tmp_path / f"h{regenotype}.vcf")
        st = H.harmonize_file(a, out, union, 0, 0.0, 0.2, regenotype=regenotype)
        assert st["alts_added"] == 0
        (alleles, info, samples), = _first(out)
        assert alleles == ("A", "C", "T")
        assert tuple(info["AD"]) == (20, 5, 10)
        assert samples["S1"][0] == (0, 2) and samples["S1"][2] == (10, 0, 10)
        assert samples["S2"][0] == (0, 1) and samples["S2"][2] == (10, 5, 0)


def test_keep_ref_only_writes_sites_nobody_varies_at(tmp_path):
    """A callset over a list of positions answers at every position: a site that is
    reference in every file is written through as REF>. instead of dropped."""
    a = _write_r(tmp_path / "a.vcf",
                 "chr1\t100\t.\tA\tT\t.\t.\tAD=25,30;DP=55\tGT:AD:ADF\t1/1:0,30:0,16\t0/0:25,0:12,0\n"
                 "chr1\t200\t.\tG\t.\t.\t.\tAD=40;DP=40\tGT:AD:ADF\t0/0:20:10\t0/0:20:10\n")
    b = _write_r(tmp_path / "b.vcf",
                 "chr1\t100\t.\tA\t.\t.\t.\tAD=50;DP=50\tGT:AD:ADF\t0/0:25:12\t0/0:25:12\n"
                 "chr1\t200\t.\tG\t.\t.\t.\tAD=40;DP=40\tGT:AD:ADF\t0/0:20:10\t0/0:20:10\n",
                 samples="S3\tS4")
    union, _d, _a, st = H.accumulate_union([a, b], 0, 0.0, 0.2)
    assert set(union) == {("chr1", 100)} and st["union_dropped"] == 1

    union, _d, _a, st = H.accumulate_union([a, b], 0, 0.0, 0.2, keep_ref_only=True)
    assert union[("chr1", 200)] == ["G"] and st["union_dropped"] == 0
    assert st["union_with_alts"] == 1
    r = H.harmonize_file(b, str(tmp_path / "hb.vcf"), union, 0, 0.0, 0.2, regenotype=False)
    assert r["written"] == 2 and r["dropped_ref_only"] == 0 and r["alts_added"] == 1
    recs = _first(str(tmp_path / "hb.vcf"))
    assert recs[0][0] == ("A", "T") and recs[0][2]["S3"] == ((0, 0), False, (25, 0), (12, 0))
    assert tuple(recs[0][1]["AD"]) == (50, 0)
    assert recs[1][0][0] == "G" and not [x for x in recs[1][0][1:] if x != "."]
    assert recs[1][2]["S3"][2] == (20,) and tuple(recs[1][1]["AD"]) == (40,)


def test_stale_format_fields_can_keep_the_padded_counts(tmp_path):
    a = _write_r(tmp_path / "a.vcf", "")
    assert H.stale_format_fields(a) == ["ADF", "PL"]
    assert H.stale_format_fields(a, keep=("GT", "AD", "ADF", "ADR")) == ["PL"]


def test_dropping_an_unsupported_allele_reshapes_every_per_allele_field(tmp_path):
    """G>T,A with A unsupported. Cleaning drops A: INFO/AD and FORMAT/ADF lose the slot
    with it, S1 (0/1 for T) is re-indexed, and only a sample that had called the dropped
    allele is re-called from AD -- when regenotype=False."""
    f = _write_r(tmp_path / "m.vcf",
                 "chr1\t100\t.\tG\tT,A\t.\t.\tAD=10,14,0;AC=2,1;DP=24\tGT:AD:ADF\t"
                 "0/1:10,5,0:5,3,0\t1/2:0,9,0:0,5,0\n")
    union, _d, _a, _st = H.accumulate_union([f], 3, 0.005, 0.2)
    assert union[("chr1", 100)] == ["G", "T"]
    H.harmonize_file(f, str(tmp_path / "out.vcf"), union, 3, 0.005, 0.2,
                     regenotype=False, stale_info=("AC",))
    (alleles, info, samples), = _first(str(tmp_path / "out.vcf"))
    assert alleles == ("G", "T")
    assert tuple(info["AD"]) == (10, 14) and "AC" not in info and info["DP"] == 24
    assert samples["S1"] == ((0, 1), False, (10, 5), (5, 3))     # kept, re-indexed
    assert samples["S2"] == ((1, 1), False, (0, 9), (0, 5))      # was 1/2, A gone: re-called


def _unsorted_vcf(tmp_path):
    """Two records at one position with another position wedged between them.

    Pass 1 collapses duplicates over the whole file, pass 2 only over *adjacent* records,
    so this is the input on which the two passes disagree about which record wins. Real
    bcftools output is coordinate-sorted and never looks like this; a hand-edited or
    concatenated file can.
    """
    return _write_vcf(tmp_path / "unsorted.vcf",
        "chr1\t384\t.\tA\tT,G\t.\t.\t.\tGT:AD\t1/1:0,26,0\t2/2:0,0,24\n"
        "chr1\t500\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"
        "chr1\t384\t.\tA\tC\t.\t.\t.\tGT:AD\t1/1:0,26\t0/0:9,0\n")


def test_a_non_adjacent_duplicate_position_is_refused_not_crashed_on(tmp_path):
    """Pass 2 must not emit a record whose alleles pass 1 left out of the union.

    Before this check the second `chr1:384` record reached `_remap_gt`, whose
    `old_to_new[g]` raised a bare `KeyError: 1` — no position, no file, no clue what was
    wrong. With `regenotype=True` there was no crash at all and the sample was silently
    re-called from an all-zero relaid AD, so its reads simply vanished.
    """
    f = _unsorted_vcf(tmp_path)
    union, _dups, _amb, _st = H.accumulate_union([f], min_ad=3, min_af=0.005,
                                                 het_min_af=0.2)
    # pass 1 preferred the 2-ALT record, so 'C' is not in the union
    assert "C" not in union[("chr1", 384)]

    out = str(tmp_path / "out.vcf")
    with pytest.raises(SystemExit, match="coordinate-sorted"):
        H.harmonize_file(f, out, union, 3, 0.005, 0.2, regenotype=False)


def test_the_same_refusal_applies_on_the_regenotype_path(tmp_path):
    """The silent variant is the more dangerous one, so it must be refused too."""
    f = _unsorted_vcf(tmp_path)
    union, _dups, _amb, _st = H.accumulate_union([f], min_ad=3, min_af=0.005,
                                                 het_min_af=0.2)
    out = str(tmp_path / "out2.vcf")
    with pytest.raises(SystemExit, match="coordinate-sorted"):
        H.harmonize_file(f, out, union, 3, 0.005, 0.2, regenotype=True)


def test_an_unmappable_allele_becomes_a_missing_call_rather_than_a_keyerror():
    """`_remap_gt`'s own guard, independent of how the record got there."""
    class _S(dict):
        phased = False
        def get(self, k, d=None):
            return dict.get(self, k, d)

    s = _S(GT=(1, 1))
    dropped = H._remap_gt(s, {0: 0, 2: 1})
    assert s["GT"] == (None, None)
    assert dropped == 2


# --- records at one position with REFs of different length ----------------------------
#
# `A > T` and `ATT > A` both sit at POS 500. With indels kept, the union used to be
# `[first file's REF] + sorted(every file's ALTs)`, which with the SNP file first came out
# as `REF=A ALT=A,T` -- an ALT equal to the reference, and every carrier of the deletion
# re-labelled `1/1` of it. Nothing errored. The alleles are only comparable once written
# against one REF, and the longer one is it: `A > T` is `ATT > TTT`, the same variant.


def test_a_snp_and_a_deletion_at_one_position_share_the_longer_ref(tmp_path):
    a = _write_vcf(tmp_path / "a.vcf",
                   "chr1\t500\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,40\t0/0:40,0\n")
    b = _write_vcf(tmp_path / "b.vcf",
                   "chr1\t500\t.\tATT\tA\t.\t.\t.\tGT:AD\t1/1:0,40\t0/0:40,0\n")
    union, _d, _a, _st = H.accumulate_union([a, b], min_ad=0, min_af=0.0, het_min_af=0.2,
                                            drop_indels=False)
    assert union[("chr1", 500)] == ["ATT", "A", "TTT"]
    assert "A" != union[("chr1", 500)][0]          # no ALT equal to REF


def test_the_carriers_end_up_on_the_right_allele_in_both_files(tmp_path):
    a = _write_vcf(tmp_path / "a.vcf",
                   "chr1\t500\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,40\t0/0:40,0\n")
    b = _write_vcf(tmp_path / "b.vcf",
                   "chr1\t500\t.\tATT\tA\t.\t.\t.\tGT:AD\t1/1:0,40\t0/0:40,0\n")
    union, _d, _a, _st = H.accumulate_union([a, b], min_ad=0, min_af=0.0, het_min_af=0.2,
                                            drop_indels=False)
    oa, ob = str(tmp_path / "ha.vcf"), str(tmp_path / "hb.vcf")
    H.harmonize_file(a, oa, union, 0, 0.0, 0.2, drop_indels=False)
    H.harmonize_file(b, ob, union, 0, 0.0, 0.2, drop_indels=False)
    (_c, _p, alleles_a, gts_a), = _records(oa)
    (_c, _p, alleles_b, gts_b), = _records(ob)
    assert alleles_a == alleles_b == ("ATT", "A", "TTT")
    assert gts_a["S1"] == (2, 2)      # the SNP carrier: TTT is index 2
    assert gts_b["S1"] == (1, 1)      # the deletion carrier: A is index 1
    assert gts_a["S2"] == gts_b["S2"] == (0, 0)


def test_refs_that_are_not_prefixes_of_each_other_are_refused(tmp_path):
    # `AT` and `AG` at one position are two different reference sequences, not two lengths
    # of one; no re-expression can reconcile them and guessing would be worse than stopping
    a = _write_vcf(tmp_path / "a.vcf",
                   "chr1\t500\t.\tAT\tA\t.\t.\t.\tGT:AD\t1/1:0,40\t0/0:40,0\n")
    b = _write_vcf(tmp_path / "b.vcf",
                   "chr1\t500\t.\tAG\tA\t.\t.\t.\tGT:AD\t1/1:0,40\t0/0:40,0\n")
    with pytest.raises(SystemExit, match="not prefixes"):
        H.accumulate_union([a, b], min_ad=0, min_af=0.0, het_min_af=0.2, drop_indels=False)


def test_a_biallelic_snp_union_is_untouched_by_the_padding_path(tmp_path):
    a = _write_vcf(tmp_path / "a.vcf",
                   "chr1\t500\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,40\t0/0:40,0\n")
    b = _write_vcf(tmp_path / "b.vcf",
                   "chr1\t500\t.\tA\tG\t.\t.\t.\tGT:AD\t1/1:0,40\t0/0:40,0\n")
    union, _d, _a, _st = H.accumulate_union([a, b], min_ad=0, min_af=0.0, het_min_af=0.2)
    assert union[("chr1", 500)] == ["A", "G", "T"]


def test_a_record_reduced_to_ref_only_does_not_invent_reference_calls(tmp_path):
    """When every ALT is cleaned away, the samples WITH reference reads become 0/0. A sample
    with no reads at all -- AD=0,0, or no AD -- has no evidence for any call, and used to be
    stamped 0/0 with the rest: the "total 0 -> missing" rule the module applies everywhere
    else, broken on this one path."""
    f = _write_vcf(tmp_path / "a.vcf",
                   "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:30,1\t./.:0,0\n")
    with pysam.VariantFile(f) as vf:
        rec = next(vf)
        H.clean_record(rec, min_ad=3, min_af=0.05, het_min_af=0.2)
        assert rec.alleles == ("A", ".")                     # the one-read ALT is gone
        assert tuple(rec.samples["S1"]["GT"]) == (0, 0)      # 30 reference reads: 0/0
        assert tuple(rec.samples["S2"]["GT"]) == (None, None)  # no reads: still missing
