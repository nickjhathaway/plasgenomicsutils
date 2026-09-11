"""WSMAF at a site with more than one alternate.

`wsmaf` is meant to be the within-sample frequency of the **population's minor allele**, so
that a sample's value can be read against the population frequency the site shows. The old
form computed the *sum of alternates* -- `1 - p_ref` -- which is that allele only when there
is one of them. At a triallelic site it is the pooled non-reference fraction, and the
`plaf > 0.5` flip that orients it presupposes a two-allele partition.

`minor_frac` was always right and is untouched: the runner-up allele's share, which is what
bounds it by 0.5 and makes it comparable to a strain proportion.
"""

import numpy as np
import pytest

pytest.importorskip("cyvcf2")

from plasgenomicsutils.lib.wsaf import wsaf_profile

_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
)


def _vcf(tmp_path, samples, rows):
    hdr = _HDR + ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(samples) + "\n")
    def cell(ad):
        # a genotype is required by the FORMAT string; wsaf reads AD, not GT, so the call
        # itself only has to be well formed
        parts = [int(x) for x in ad.split(",")]
        k = max(range(len(parts)), key=lambda i: parts[i])
        return f"{k}/{k}:{ad}"

    body = "".join(f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT:AD\t"
                   + "\t".join(cell(c) for c in cells) + "\n"
                   for pos, alt, cells in rows)
    p = tmp_path / "in.vcf"
    p.write_text(hdr + body)
    return str(p)


def _sites(path):
    import csv
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_wsmaf_is_the_named_minor_alleles_frequency_not_the_sum_of_alternates(tmp_path):
    """Two alternates, one common and one rare, in a sample carrying both.

    Reference 20 reads, C 60, G 20. The population's minor allele here is the reference at
    20%, so the sample's wsmaf should be its reference fraction, 0.2 -- not `1 - p_ref`
    (0.8) and not the pooled non-reference share.
    """
    samples = ["s1", "s2", "s3", "s4"]
    rows = [(1000, "C,G", ["20,60,20", "0,100,0", "0,100,0", "0,100,0"])]
    out = str(tmp_path / "sites.tsv")
    wsaf_profile(_vcf(tmp_path, samples, rows), min_minor=0.05, min_minor_reads=2,
                 sites_out=out)
    rec = [r for r in _sites(out) if r["sample"] == "s1"]
    assert len(rec) == 1
    assert float(rec[0]["wsmaf"]) == pytest.approx(0.2, abs=1e-6)


def test_the_minor_fraction_is_still_the_runner_up_share(tmp_path):
    """Untouched: three alleles at a third each read as 0.33, not 0.67."""
    samples = ["s1", "s2"]
    rows = [(1000, "C,G", ["30,30,30", "0,90,0"])]
    out = str(tmp_path / "sites.tsv")
    wsaf_profile(_vcf(tmp_path, samples, rows), min_minor=0.05, min_minor_reads=2,
                 sites_out=out)
    rec = [r for r in _sites(out) if r["sample"] == "s1"][0]
    assert float(rec["minor_frac"]) == pytest.approx(1 / 3, abs=1e-6)


def test_a_biallelic_site_is_unchanged(tmp_path):
    """Where the old reading was right, the new one must give the same number."""
    samples = ["s1", "s2", "s3", "s4"]
    # reference is the minor allele overall (25 reads of 100 pooled), so wsmaf is the
    # sample's reference fraction
    rows = [(1000, "C", ["70,30", "0,100", "0,100", "0,100"])]
    out = str(tmp_path / "sites.tsv")
    wsaf_profile(_vcf(tmp_path, samples, rows), min_minor=0.05, min_minor_reads=2,
                 sites_out=out)
    rec = [r for r in _sites(out) if r["sample"] == "s1"][0]
    assert float(rec["wsmaf"]) == pytest.approx(0.7, abs=1e-6)


def test_the_site_table_names_which_allele_wsmaf_is_about(tmp_path):
    """A number that means "this allele's frequency" has to say which allele."""
    samples = ["s1", "s2", "s3", "s4"]
    rows = [(1000, "C,G", ["20,60,20", "0,100,0", "0,100,0", "0,100,0"])]
    out = str(tmp_path / "sites.tsv")
    wsaf_profile(_vcf(tmp_path, samples, rows), min_minor=0.05, min_minor_reads=2,
                 sites_out=out)
    rec = [r for r in _sites(out) if r["sample"] == "s1"][0]
    assert "wsmaf_allele" in rec
    assert rec["wsmaf_allele"] == "A"       # the reference is the population minor here
