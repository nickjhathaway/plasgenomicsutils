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


# ---- SNP-ness and allele count are separate questions ------------------------------

def _types_vcf(tmp_path):
    """One record of each shape the two tests have to tell apart."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"]
    rows = [(100, "A", "T", "0/1", "1/1"),          # biallelic SNP
            (200, "A", "T,G", "1/2", "1/1"),        # multiallelic SNP
            (300, "A", "T,ATT", "1/2", "1/1"),      # mixed SNP + indel
            (400, "ATT", "A", "0/1", "1/1"),        # indel
            (500, "A", ".", "0/0", "0/0"),          # no ALT
            (600, "AT", "GC", "0/1", "1/1")]        # MNP
    for pos, ref, alt, g1, g2 in rows:
        hdr.append(f"c1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{g1}\t{g2}")
    p = tmp_path / "types.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _kept(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\n", path], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    return [int(x) for x in out.stdout.split()]


@pytest.mark.parametrize("kw,expected", [
    ({}, [100]),                                        # biallelic SNPs
    ({"biallelic": False}, [100, 200]),                 # SNPs, multiallelic allowed
    ({"snps_only": False}, [100, 400, 600]),            # biallelic anything
])
def test_the_two_tests_are_independent(tmp_path, kw, expected):
    """`mnp_handling` is pinned here so this stays a test of the other two knobs."""
    out = str(tmp_path / "o.bcf")
    F.biallelic_snp_filter(_types_vcf(tmp_path), out, trim=False,
                           mnp_handling="remove", **kw)
    assert _kept(out) == expected


def test_a_mixed_snp_indel_record_is_not_a_snp(tmp_path):
    """`bcftools view -v snps` keeps a record if ANY allele is a SNP, so `A>T,ATT` passes it.
    Selecting by exclusion is what makes 'SNPs only' mean every allele."""
    out = str(tmp_path / "o.bcf")
    F.biallelic_snp_filter(_types_vcf(tmp_path), out, trim=False, biallelic=False,
                           mnp_handling="remove")
    assert 300 not in _kept(out)


def test_a_record_left_with_no_alt_is_not_kept_by_the_type_test(tmp_path):
    """With the biallelic test off nothing else drops it, so `ref` has to be excluded too."""
    out = str(tmp_path / "o.bcf")
    F.biallelic_snp_filter(_types_vcf(tmp_path), out, trim=False, biallelic=False,
                           mnp_handling="remove")
    assert 500 not in _kept(out)


def test_trimming_still_rescues_a_site_that_becomes_biallelic(tmp_path):
    """The reason the step runs after re-genotyping, unchanged by the split."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
           "c1\t100\t.\tA\tT,G\t.\t.\t.\tGT\t0/1\t1/1"]     # nothing carries the G
    v = tmp_path / "trim.vcf"
    v.write_text("\n".join(hdr) + "\n")
    trimmed, plain = str(tmp_path / "a.bcf"), str(tmp_path / "b.bcf")
    F.biallelic_snp_filter(str(v), trimmed, trim=True)
    F.biallelic_snp_filter(str(v), plain, trim=False)
    assert _kept(trimmed) == [100] and _kept(plain) == []


def test_turning_everything_off_is_refused_rather_than_copying_the_input(tmp_path):
    with pytest.raises(SystemExit, match="drop it from the chain"):
        F.biallelic_snp_filter(_types_vcf(tmp_path), str(tmp_path / "o.bcf"),
                               trim=False, snps_only=False, biallelic=False,
                               mnp_handling="keep")


def test_with_both_tests_off_it_trims_and_says_so(tmp_path, capsys):
    """Legal, since trimming is real work -- but it must not look like a filter that ran."""
    out = str(tmp_path / "o.bcf")
    F.biallelic_snp_filter(_types_vcf(tmp_path), out, trim=True, snps_only=False,
                           biallelic=False, mnp_handling="keep")
    assert "only trims unused ALT alleles" in capsys.readouterr().out
    assert set(_kept(out)) == {100, 200, 300, 400, 500, 600}


@pytest.mark.parametrize("snps_only,biallelic", [(True, True), (True, False), (False, True)])
def test_the_type_test_excludes_a_no_alt_record_by_allele_count(snps_only, biallelic):
    """Not by trusting `-V ref`, whose meaning depends on the bcftools version.

    A record with ALT="." is type `ref` to bcftools 1.24 but to 1.19 is no type at all --
    `-v ref` selects nothing there, so `-V ref` excludes nothing and a non-variant record
    survives a SNPs-only filter. `--min-alleles 2` counts alleles instead, which every
    version agrees on, so it has to be present whenever either test is on.
    """
    args = F._snp_select_args(snps_only, biallelic)
    assert "-m2" in args.split()
    assert ("-M2" in args.split()) == biallelic
    assert (f"-V {F.NON_SNP_TYPES}" in args) == snps_only


# --- mnp_handling ----------------------------------------------------------------------

def _mnp_vcf(tmp_path):
    """A plain SNP, a true MNP, and a SNP written with padding."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
           "c1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1",
           "c1\t200\t.\tAT\tGC\t.\t.\t.\tGT\t0/1\t1/1",       # a real MNP
           "c1\t300\t.\tTTATA\tCTATA\t.\t.\t.\tGT\t0/1\t1/1"]  # one base differs
    p = tmp_path / "mnp.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _rec(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\t%REF\t%ALT\n", path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return [tuple(line.split("\t")) for line in out.stdout.splitlines()]


def test_split_breaks_an_mnp_into_snps_and_minimises_a_padded_one(tmp_path):
    """The default. `TTATA>CTATA` is one substitution written long, not an MNP -- splitting
    rewrites it minimally, which is what a downstream tool assuming one base per record
    needs."""
    out = str(tmp_path / "o.bcf")
    F.biallelic_snp_filter(_mnp_vcf(tmp_path), out, trim=False)
    assert _rec(out) == [("100", "A", "T"), ("200", "A", "G"), ("201", "T", "C"),
                         ("300", "T", "C")]


def test_remove_drops_an_mnp_and_keep_leaves_it(tmp_path):
    dropped, kept = str(tmp_path / "d.bcf"), str(tmp_path / "k.bcf")
    F.biallelic_snp_filter(_mnp_vcf(tmp_path), dropped, trim=False, mnp_handling="remove")
    F.biallelic_snp_filter(_mnp_vcf(tmp_path), kept, trim=False, mnp_handling="keep")
    # the padded record is a SNP to bcftools either way -- only the real MNP moves
    assert ("200", "AT", "GC") not in _rec(dropped)
    assert ("200", "AT", "GC") in _rec(kept)
    assert ("300", "TTATA", "CTATA") in _rec(kept)


def test_split_leaves_multiallelic_sites_whole(tmp_path):
    """`norm -a` decomposes multiallelic sites as readily as MNPs, one record per ALT with
    `*` for the others. Splitting by allele count first keeps it away from the sites that
    `biallelic=False` exists to preserve."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
           "c1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1",
           "c1\t200\t.\tA\tT,G\t.\t.\t.\tGT\t1/2\t1/1",     # multiallelic SNP
           "c1\t300\t.\tAT\tGC\t.\t.\t.\tGT\t0/1\t1/1"]     # MNP
    v = tmp_path / "both.vcf"
    v.write_text("\n".join(hdr) + "\n")
    out = str(tmp_path / "o.bcf")
    F.biallelic_snp_filter(str(v), out, trim=False, biallelic=False, mnp_handling="split")
    # the MNP is split, the multiallelic site is untouched, and the order is restored
    assert _rec(out) == [("100", "A", "T"), ("200", "A", "T,G"),
                         ("300", "A", "G"), ("301", "T", "C")]


def test_a_record_whose_only_alt_is_a_spanning_deletion_is_dropped(tmp_path):
    """`T > *` says an upstream deletion covers this position and nothing else: two alleles
    by count, no variant to call. `A > *,T` is a real SNP that merely sits under one."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
           "c1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1",
           "c1\t200\t.\tT\t*\t.\t.\t.\tGT\t0/1\t1/1",
           "c1\t300\t.\tA\t*,T\t.\t.\t.\tGT\t1/2\t1/1"]
    v = tmp_path / "star.vcf"
    v.write_text("\n".join(hdr) + "\n")
    strict, loose = str(tmp_path / "a.bcf"), str(tmp_path / "b.bcf")
    F.biallelic_snp_filter(str(v), strict, trim=False)
    F.biallelic_snp_filter(str(v), loose, trim=False, biallelic=False,
                           mnp_handling="keep")
    assert _kept(strict) == [100]                    # the multiallelic goes to `biallelic`
    assert _kept(loose) == [100, 300]                # ...and comes back when it is off
    assert 200 not in _kept(loose)


def test_an_unknown_mnp_handling_is_named(tmp_path):
    with pytest.raises(SystemExit, match="mnp_handling must be one of"):
        F.biallelic_snp_filter(_mnp_vcf(tmp_path), str(tmp_path / "o.bcf"),
                               mnp_handling="atomise")


def test_split_breaks_up_a_multiallelic_mnp_too(tmp_path):
    """Splitting the site apart first is what makes this reachable: atomising a multiallelic
    record directly would decompose the site instead of the substitution."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=10000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
           "c1\t300\t.\tAT\tGC,GG\t.\t.\t.\tGT\t1/2\t1/1"]
    v = tmp_path / "mm.vcf"
    v.write_text("\n".join(hdr) + "\n")
    out = str(tmp_path / "o.bcf")
    F.biallelic_snp_filter(str(v), out, trim=False, biallelic=False, mnp_handling="split")
    # both ALTs carry G at the first base, so that position collapses; the second stays
    # multiallelic because that is where the two alleles actually differ
    assert _rec(out) == [("300", "A", "G"), ("301", "T", "C,G")]
