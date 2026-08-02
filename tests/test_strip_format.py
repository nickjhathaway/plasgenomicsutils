"""strip_stale_format: null genotype-linked FORMAT fields (e.g. PL) whose length no
longer matches the genotypes, and the filter_ad_regenotype self-clean that prevents the
inconsistency in the first place."""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib.strip_format import strip_stale_format

HAVE_BCFTOOLS = shutil.which("bcftools") is not None
needs = pytest.mark.skipif(not HAVE_BCFTOOLS, reason="bcftools not on PATH")
pysam = pytest.importorskip("pysam")

# 100: triallelic A>T,G, G unused (trims to biallelic); S2 has a 7-value PL (a diploid GT
#      forced over higher-ploidy calls) -> inconsistent. 200: an ordinary biallelic SNP.
_VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="phred GLs">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
    "chr1\t100\t.\tA\tT,G\t.\t.\t.\tGT:PL\t0/1:0,10,20,30,40,50\t0/0:0,10,20,30,40,50,60\n"
    "chr1\t200\t.\tC\tG\t.\t.\t.\tGT:PL\t0/1:0,10,20\t1/1:30,10,0\n"
)


def _bgzip(tmp_path):
    p = tmp_path / "in.vcf"
    p.write_text(_VCF)
    subprocess.run(["bgzip", "-f", str(p)], check=True)
    return str(p) + ".gz"


def _pl(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS[ %PL]\n", path],
                         stdout=subprocess.PIPE, text=True).stdout
    return {r.split()[0]: r.split()[1:] for r in out.splitlines()}


def _trims_ok(path):
    r = subprocess.run("bcftools view --trim-alt-alleles '%s' -Ou | "
                       "bcftools view -m2 -M2 -v snps - -Ob -o /dev/null" % path,
                       shell=True, stderr=subprocess.DEVNULL)
    return r.returncode == 0


@needs
def test_mismatch_mode_is_surgical(tmp_path):
    inp = _bgzip(tmp_path)
    out = str(tmp_path / "out.bcf")
    n = strip_stale_format(inp, out, fields=("PL",), mode="mismatch")
    assert n == 1                                   # only the one inconsistent record touched
    pl = _pl(out)
    assert pl["200"] == ["0,10,20", "30,10,0"]      # valid record's PL preserved
    assert set(pl["100"]) == {".,.,.,.,.,."}        # inconsistent record nulled (G=6, diploid)
    assert _trims_ok(out)                           # and bcftools can now trim


@needs
def test_always_mode_drops_the_field(tmp_path):
    inp = _bgzip(tmp_path)
    out = str(tmp_path / "out.bcf")
    strip_stale_format(inp, out, fields=("PL",), mode="always")
    hdr = subprocess.run(["bcftools", "view", "-h", out], stdout=subprocess.PIPE, text=True).stdout
    assert "ID=PL" not in hdr                        # field removed entirely
    assert _trims_ok(out)


@needs
def test_clean_file_is_untouched(tmp_path):
    # a file with no PL: mismatch mode is a no-op passthrough and stays trimmable
    p = tmp_path / "clean.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/1\n")
    subprocess.run(["bgzip", "-f", str(p)], check=True)
    out = str(tmp_path / "out.bcf")
    assert strip_stale_format(str(p) + ".gz", out, mode="mismatch") == 0
    assert _positions_ok(out)


def _positions_ok(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\n", path],
                         stdout=subprocess.PIPE, text=True).stdout
    return out.split() == ["10"]


def _regen_input(tmp_path, body, header_extra=""):
    inp = tmp_path / "in.vcf"
    inp.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        '##FORMAT=<ID=ADS,Number=1,Type=Integer,Description="sum">\n'
        + header_extra +
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n" + body)
    subprocess.run(["bgzip", "-f", str(inp)], check=True)
    return str(inp) + ".gz"


@needs
def test_ploidy_reduce_to_haploid_warns_and_trims(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.regenotype import filter_ad_regenotype
    inp = _regen_input(
        tmp_path,
        "chr1\t100\t.\tA\tT,G\t.\t.\t.\tGT:AD:ADS:PL"
        "\t0/1:5,5,0:10:0,10,20,30,40,50\t0/0:9,0,0:9:0,10,20,30,40,50,60\n",
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="GL">\n')
    out = str(tmp_path / "out.bcf")
    filter_ad_regenotype(inp, out, ploidy=1)
    gt = subprocess.run(["bcftools", "query", "-f", "[%GT ]\n", out],
                        stdout=subprocess.PIPE, text=True).stdout.split()
    assert all("/" not in g and "|" not in g for g in gt)   # haploid GTs
    # PL trimmed to haploid Number=G = n_alleles (3 for this triallelic record), all missing
    pl = _pl(out)["100"]
    assert all(v.count(",") + 1 == 3 for v in pl)


@needs
def test_ploidy_promotion_errors(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.regenotype import filter_ad_regenotype
    inp = _regen_input(tmp_path, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD:ADS\t0:9,0:9\t1:0,9:9\n")
    with pytest.raises(SystemExit):
        filter_ad_regenotype(inp, str(tmp_path / "out.bcf"), ploidy=2)  # 2 > haploid input


@needs
def test_default_keeps_diploid_coding_on_haploid_input(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.regenotype import filter_ad_regenotype
    inp = _regen_input(tmp_path, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD:ADS\t0:9,0:9\t1:0,9:9\n")
    out = str(tmp_path / "out.bcf")
    filter_ad_regenotype(inp, out)   # no ploidy -> conventional diploid coding, no error
    gt = subprocess.run(["bcftools", "query", "-f", "[%GT ]\n", out],
                        stdout=subprocess.PIPE, text=True).stdout.split()
    assert all(len(g.replace("|", "/").split("/")) == 2 for g in gt)   # diploid


@needs
def test_regenotype_self_cleans_stale_likelihoods(tmp_path):
    pytest.importorskip("cyvcf2")
    from plasgenomicsutils.lib.regenotype import filter_ad_regenotype
    # AD/ADS present so the record is re-genotyped (to diploid); the 7-wide PL on S2 must be
    # blanked to a diploid-consistent length so downstream trimming works.
    inp = tmp_path / "in.vcf"
    inp.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        '##FORMAT=<ID=ADS,Number=1,Type=Integer,Description="sum">\n'
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="GL">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "chr1\t100\t.\tA\tT,G\t.\t.\t.\tGT:AD:ADS:PL"
        "\t0/1:5,5,0:10:0,10,20,30,40,50\t0/0:9,0,0:9:0,10,20,30,40,50,60\n")
    subprocess.run(["bgzip", "-f", str(inp)], check=True)
    out = str(tmp_path / "out.bcf")
    filter_ad_regenotype(str(inp) + ".gz", out)
    # PL blanked (all missing) and the record now trims cleanly
    pl = _pl(out)["100"]
    assert all(set(v) <= {".", ","} for v in pl)
    assert _trims_ok(out)
