"""Tests for the region filters and the biallelic-SNP (trim) step.

These shell out to bcftools/bedtools; they skip cleanly where those aren't installed.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import vcf_filters as F
from plasgenomicsutils.lib import assets

HAVE_BCFTOOLS = shutil.which("bcftools") is not None
HAVE_BEDTOOLS = shutil.which("bedtools") is not None
needs_bcftools = pytest.mark.skipif(not HAVE_BCFTOOLS, reason="bcftools not on PATH")
needs_both = pytest.mark.skipif(not (HAVE_BCFTOOLS and HAVE_BEDTOOLS),
                                reason="bcftools+bedtools not on PATH")

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)


def _vcf(path, body):
    path.write_text(HEADER + body)
    return str(path)


def _count(path):
    out = subprocess.run(["bcftools", "view", "-H", path], stdout=subprocess.PIPE, text=True)
    return len([x for x in out.stdout.splitlines() if x])


def _positions(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\n", path], stdout=subprocess.PIPE, text=True)
    return sorted(int(x) for x in out.stdout.split())


@needs_bcftools
def test_biallelic_snp_filter_trims_artifact_allele(tmp_path):
    # 100: clean biallelic SNP (kept)
    # 200: "multiallelic" A>T,G but G is in NO genotype -> trimmed -> biallelic -> kept
    # 300: genuinely multiallelic (both T and G used) -> stays multiallelic -> dropped
    # 400: indel -> dropped (not a SNP)
    inp = _vcf(tmp_path / "in.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t1/1:0,9\n"
        "chr1\t200\t.\tA\tT,G\t.\t.\t.\tGT:AD\t0/1:5,5,0\t1/1:0,9,0\n"
        "chr1\t300\t.\tA\tT,G\t.\t.\t.\tGT:AD\t0/1:5,5,0\t2/2:0,0,9\n"
        "chr1\t400\t.\tAT\tA\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n")
    out = str(tmp_path / "out.vcf")
    F.biallelic_snp_filter(inp, out, trim=True)
    assert _positions(out) == [100, 200]   # 200 rescued by trim; 300 (multi) & 400 (indel) dropped

    # without trim, 200 stays multiallelic and is dropped
    out2 = str(tmp_path / "out2.vcf")
    F.biallelic_snp_filter(inp, out2, trim=False)
    assert _positions(out2) == [100]


_PL_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="phred-scaled GLs">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)


@needs_bcftools
def test_biallelic_trim_survives_inconsistent_genotype_likelihoods(tmp_path):
    # 100: triallelic A>T,G with G in no genotype (trimmed to biallelic), but S2 carries a
    # PL of 7 values (a mixed-ploidy / stale-likelihood artifact). Raw
    # `bcftools view --trim-alt-alleles` aborts on the Number=G length mismatch
    # ("expected 3, found 7"); the step must drop the stale likelihoods first and still
    # rescue the biallelic SNP. 200 is an ordinary biallelic SNP.
    inp = tmp_path / "in.vcf"
    inp.write_text(_PL_HEADER +
        "chr1\t100\t.\tA\tT,G\t.\t.\t.\tGT:PL\t0/1:0,10,20,30,40,50\t0/0:0,10,20,30,40,50,60\n"
        "chr1\t200\t.\tC\tG\t.\t.\t.\tGT:PL\t0/1:0,10,20\t1/1:30,10,0\n")
    out = str(tmp_path / "out.bcf")
    F.biallelic_snp_filter(str(inp), out, trim=True)   # must not raise
    assert _positions(out) == [100, 200]
    # surgical: the field is kept; only the inconsistent record (100) is nulled, while the
    # consistent record (200) keeps its real PL
    pl = subprocess.run(["bcftools", "query", "-f", "%POS[ %PL]\n", out],
                        stdout=subprocess.PIPE, text=True).stdout
    rows = {r.split()[0]: r.split()[1:] for r in pl.splitlines()}
    assert rows["200"] == ["0,10,20", "30,10,0"]       # valid likelihoods preserved
    assert set(rows["100"]) == {".,.,."}               # inconsistent record nulled (trimmed to G=3)


_GRP_HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n"
)


@needs_bcftools
def test_maf_filter_grouped_keeps_union_and_preserves_genotypes(tmp_path, bgzip):
    # groups A = {s1,s2}, B = {s3,s4}
    # 100: monomorphic in A (0/0) but polymorphic in B -> kept, A's genotypes preserved
    # 200: monomorphic everywhere -> dropped
    # 300: polymorphic only in A -> kept
    p = tmp_path / "in.vcf"
    p.write_text(_GRP_HEADER +
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\t0/1\t0/1\n"
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0\n"
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\t0/0\t0/0\n")
    inp = bgzip(p)
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tgroup\ns1\tA\ns2\tA\ns3\tB\ns4\tB\n")

    out = str(tmp_path / "out.bcf")
    F.maf_filter(inp, out, maf_min=0.02, meta=str(meta), group_col="group")
    assert _positions(out) == [100, 300]                       # union of per-group passers
    samples = subprocess.run(["bcftools", "query", "-l", out],
                             stdout=subprocess.PIPE, text=True).stdout.split()
    assert samples == ["s1", "s2", "s3", "s4"]                 # every sample retained
    # site 100 kept even though it is monomorphic in group A, and A's 0/0 calls are intact
    gt = subprocess.run(["bcftools", "query", "-f", "%POS[ %GT]\n", out],
                        stdout=subprocess.PIPE, text=True).stdout.splitlines()
    row100 = next(r for r in gt if r.startswith("100"))
    assert row100.split()[1:] == ["0/0", "0/0", "0/1", "0/1"]  # not emptied to "."


@needs_bcftools
def test_maf_filter_grouped_preserves_subthreshold_carrier_genotypes(tmp_path, bgzip):
    # 60 samples in group A + 4 in group B. a01 carries the alt but its frequency is < 2%
    # WITHIN group A; the site is rescued by group B, and a01's real 0/1 / 1/1 calls survive.
    A = [f"a{i:02d}" for i in range(1, 61)]
    B = [f"b{i}" for i in range(1, 5)]
    samples = A + B
    hdr = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")

    def row(pos, gts):
        return f"chr1\t{pos}\t.\tA\tT\t.\t.\t.\tGT\t" + "\t".join(gts.get(s, "0/0") for s in samples) + "\n"

    body = (row(100, {"a01": "0/1", "b1": "0/1", "b2": "0/1"})    # a01 het,  MAF within A ~0.008
            + row(200, {"a01": "1/1", "b1": "0/1", "b2": "0/1"})  # a01 hom,  MAF within A ~0.017
            + row(300, {"a01": "0/1"}))                            # rare in EVERY group -> dropped
    p = tmp_path / "in.vcf"
    p.write_text(hdr + body)
    p = bgzip(p)
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tcountry\n"
                    + "".join(f"{s}\tA\n" for s in A) + "".join(f"{s}\tB\n" for s in B))
    out = str(tmp_path / "out.bcf")
    F.maf_filter(p, out, maf_min=0.02, meta=str(meta), group_col="country")
    assert _positions(out) == [100, 200]                          # 300 (rare everywhere) dropped
    gt = dict(r.split() for r in subprocess.run(
        ["bcftools", "query", "-s", "a01", "-f", "%POS [%GT]\n", out],
        stdout=subprocess.PIPE, text=True).stdout.splitlines())
    assert gt["100"] == "0/1"   # het carrier kept, though < 2% within its own country
    assert gt["200"] == "1/1"   # hom-alt carrier kept too


@needs_bcftools
def test_snp_bed_matches_vcf_snp_ids(tmp_path):
    from plasgenomicsutils.lib.vcf_io import SnpPanel
    inp = _vcf(tmp_path / "in.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t1/1:0,9\n"     # SNP
        "chr1\t250\t.\tC\tG\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"     # SNP
        "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n")   # indel -> excluded
    bed = str(tmp_path / "out.snps.bed")
    F.snp_bed(inp, bed)
    rows = [ln.split("\t") for ln in open(bed).read().splitlines()]
    assert [r[1:] for r in rows] == [["99", "100", "chr1:99"], ["249", "250", "chr1:249"]]
    # loads as a SNP panel identical to the VCF-derived ids/coords (build_ibd_matrix input)
    p = SnpPanel.from_bed(bed).df
    assert list(p.snp_id) == ["chr1:99", "chr1:249"]
    assert list(p.pos0) == [99, 249]


@needs_bcftools
def test_maf_filter_max_defaults_to_symmetric_window(tmp_path):
    # 2 diploid samples (4 alleles). AF: 100 = 0.5, 200 = 0.25, 300 = 0.75.
    inp = _vcf(tmp_path / "in.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/1:5,5\n"     # AF 0.50
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"     # AF 0.25
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t1/1:0,9\n")    # AF 0.75
    # maf_max unset -> 1 - maf_min = 0.7, so the [0.3, 0.7] window keeps only AF 0.5
    out = str(tmp_path / "out.vcf")
    F.maf_filter(inp, out, maf_min=0.3)
    assert _positions(out) == [100]
    # an explicit asymmetric window is still honoured
    out2 = str(tmp_path / "out2.vcf")
    F.maf_filter(inp, out2, maf_min=0.3, maf_max=0.8)   # now 0.75 also passes
    assert _positions(out2) == [100, 300]


@needs_bcftools
def test_singleton_add_ads_correct_on_multiallelic(tmp_path):
    # ADS must be the sum of ALL AD entries (ref + every ALT), not just ref+alt1.
    inp = _vcf(tmp_path / "in.vcf",
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:10,6\t1/1:0,9\n"        # biallelic: 16, 9
        "chr1\t200\t.\tA\tT,G\t.\t.\t.\tGT:AD\t1/2:2,7,5\t0/1:4,8,0\n")  # triallelic: 14, 12
    out = str(tmp_path / "out.vcf")
    F.singleton_add_ads(inp, out, min_samples=0)  # keep both sites
    rows = subprocess.run(["bcftools", "query", "-f", "%POS[ %ADS]\n", out],
                          stdout=subprocess.PIPE, text=True).stdout.strip().splitlines()
    ads = {r.split()[0]: r.split()[1:] for r in rows}
    assert ads["100"] == ["16", "9"]
    assert ads["200"] == ["14", "12"]   # ref+alt1+alt2, not the old ref+alt1 (would be 9 and 12)


@needs_both
def test_region_filter_keep_and_exclude(tmp_path):
    inp = _vcf(tmp_path / "in.vcf",
        "chr1\t150\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"
        "chr1\t550\t.\tA\tT\t.\t.\t.\tGT:AD\t1/1:0,9\t0/0:9,0\n"
        "chr1\t950\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n")
    bed = tmp_path / "region.bed"
    bed.write_text("chr1\t100\t200\nchr1\t900\t1000\n")  # covers 150 and 950, not 550

    keep = str(tmp_path / "keep.vcf")
    F.region_filter(inp, keep, bed=str(bed), exclude=False)
    assert _positions(keep) == [150, 950]

    drop = str(tmp_path / "drop.vcf")
    F.region_filter(inp, drop, bed=str(bed), exclude=True)
    assert _positions(drop) == [550]


def test_bundled_assets_exist():
    for name in ("pf3d7_core_regions", "pf3d7_paralog_genes", "pf3d7_tandem_repeats"):
        assert name in assets.available_assets()
        p = assets.asset_path("builtin:" + name)
        with open(p) as fh:
            assert fh.readline()  # non-empty
    assert assets.resolve_bed("/some/plain/path.bed") == "/some/plain/path.bed"
    with pytest.raises(SystemExit):
        assets.asset_path("builtin:does_not_exist")
