"""Integration test for the cyvcf2-backed filter_ad_regenotype."""

import subprocess

import pytest

pytest.importorskip("cyvcf2")
from plasgenomicsutils.lib.regenotype import filter_ad_regenotype

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    '##FORMAT=<ID=ADS,Number=1,Type=Integer,Description="summed AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
)


def _query(path, fmt):
    out = subprocess.run(["bcftools", "query", "-f", f"%POS[ {fmt}]\n", path],
                         stdout=subprocess.PIPE, text=True).stdout
    return {r.split()[0]: r.split()[1:] for r in out.splitlines()}


@pytest.mark.skipif(subprocess.run(["bash", "-c", "command -v bcftools"],
                                   stdout=subprocess.DEVNULL).returncode != 0,
                    reason="bcftools not on PATH")
def test_filter_ad_regenotype_cleans_and_regenotypes(tmp_path):
    # S1 100: alt has 1 read (< min_reads 2) -> zeroed -> 0/0, ADS 30->29
    # S2 100: balanced 15,15 -> stays het 0/1
    # S3 100: hom-alt 0,40 -> stays 1/1
    inp = tmp_path / "in.vcf"
    inp.write_text(
        HEADER +
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD:ADS\t0/1:29,1:30\t0/1:15,15:30\t1/1:0,40:40\n")
    out = str(tmp_path / "out.vcf")
    filter_ad_regenotype(str(inp), out, min_reads=2, min_freq=0.01, het_min_af=0.2)

    gt = _query(out, "%GT")
    ad = _query(out, "%AD")
    ads = _query(out, "%ADS")
    assert gt["100"] == ["0/0", "0/1", "1/1"]        # S1 artifact allele dropped -> hom-ref
    assert ad["100"] == ["29,0", "15,15", "0,40"]    # S1 alt zeroed
    assert ads["100"] == ["29", "30", "40"]          # S1 ADS recomputed from cleaned AD
