"""`filter_ad_regenotype --singletons`: apply the singleton rule to the genotypes the
re-genotyping just produced, without adding a step to the chain.

The earlier `singleton_filter_add_ads` judges carrier counts on the caller's genotypes.
Re-genotyping then moves them: calls go missing for thin depth, others are called afresh
from the cleaned AD. An allele that had two carriers there can have one here -- and for a
`*` that costs the whole record at `biallelic_snp_filter --snps-only`, taking a SNP the rest
of the cohort carries out of the panel with it.
"""

import shutil

import pytest

pysam = pytest.importorskip("pysam")
pytest.importorskip("cyvcf2")

from plasgenomicsutils.lib.regenotype import filter_ad_regenotype

_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
)
_SAMPLES = [f"s{i}" for i in range(1, 7)]


def _vcf(tmp_path, rows, name="in.vcf"):
    """rows: (pos, alt, [(gt, [ad...]), ...])"""
    hdr = _HDR + ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(_SAMPLES) + "\n")
    body = ""
    for pos, alt, cells in rows:
        body += (f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT:AD\t"
                 + "\t".join(f"{gt}:" + ",".join(str(x) for x in ad) for gt, ad in cells)
                 + "\n")
    p = tmp_path / name
    p.write_text(hdr + body)
    return str(p)


def _run(tmp_path, rows, **kw):
    out = str(tmp_path / "out.vcf")
    filter_ad_regenotype(_vcf(tmp_path, rows), out, **kw)
    with pysam.VariantFile(out) as f:
        return [(r.pos, list(r.alts or []),
                 [tuple(r.samples[s]["GT"]) for s in _SAMPLES]) for r in f]


# `A > T,*`: five samples carry T on solid depth, one carries the deletion. Exactly the
# shape that leaves the SNP panel for one sample's sake.
_SNP_WITH_LONE_STAR = [(1000, "T,*", [("1/1", [0, 40, 0])] * 5 + [("2/2", [0, 0, 40])])]


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_star_mode_blanks_the_lone_deletion_and_trims_it(tmp_path):
    (pos, alts, gts), = _run(tmp_path, _SNP_WITH_LONE_STAR, singletons="star")
    assert alts == ["T"]                        # a plain biallelic SNP now
    assert gts[:5] == [(1, 1)] * 5              # the cohort's calls are untouched
    assert gts[5] == (None, None)               # the one deleted sample is missing


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_star_is_the_default(tmp_path):
    """On by default: the step that creates the singleton `*` is the one that cleans it up,
    so a caller who does not know the problem exists does not ship a panel missing the SNP."""
    (pos, alts, gts), = _run(tmp_path, _SNP_WITH_LONE_STAR)
    assert alts == ["T"] and gts[5] == (None, None)


@pytest.mark.parametrize("off", [None, "none"])
def test_it_can_be_turned_off(tmp_path, off):
    (pos, alts, gts), = _run(tmp_path, _SNP_WITH_LONE_STAR, singletons=off)
    assert alts == ["T", "*"] and gts[5] == (2, 2)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_a_deletion_above_the_threshold_is_left_for_spanning_del_filter(tmp_path):
    rows = [(1000, "T,*", [("1/1", [0, 40, 0])] * 4 + [("2/2", [0, 0, 40])] * 2)]
    (pos, alts, gts), = _run(tmp_path, rows, singletons="star")
    assert alts == ["T", "*"] and gts[4] == (2, 2) and gts[5] == (2, 2)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_the_threshold_is_configurable(tmp_path):
    rows = [(1000, "T,*", [("1/1", [0, 40, 0])] * 4 + [("2/2", [0, 0, 40])] * 2)]
    (pos, alts, gts), = _run(tmp_path, rows, singletons="star", singleton_min_samples=2)
    assert alts == ["T"] and gts[4] == (None, None) and gts[5] == (None, None)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_star_mode_does_not_touch_a_private_real_allele(tmp_path):
    rows = [(1000, "T,G", [("1/1", [0, 40, 0])] * 5 + [("2/2", [0, 0, 40])])]
    (pos, alts, gts), = _run(tmp_path, rows, singletons="star")
    assert alts == ["T", "G"] and gts[5] == (2, 2)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_real_and_all_modes_take_the_private_base(tmp_path):
    rows = [(1000, "T,G", [("1/1", [0, 40, 0])] * 5 + [("2/2", [0, 0, 40])])]
    for mode in ("real", "all"):
        (pos, alts, gts), = _run(tmp_path, rows, singletons=mode)
        assert alts == ["T"], mode
        assert gts[5] == (None, None), mode


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_the_count_is_taken_after_re_genotyping_not_before(tmp_path):
    """The whole point. Two samples' genotypes name the `*` on the way in, but one of them
    has a single supporting read, which the AD cleaning zeroes. After re-genotyping the
    deletion is private to one sample, and only a rule applied here can see that."""
    rows = [(1000, "T,*", [("1/1", [0, 40, 0])] * 4
             + [("2/2", [0, 0, 40]), ("2/2", [0, 39, 1])])]
    (pos, alts, gts), = _run(tmp_path, rows, singletons="star",
                             min_reads=2, min_freq=0.01)
    assert gts[5] == (1, 1)                     # re-genotyped off the deletion
    assert alts == ["T"]                        # ...leaving it private, so it goes
    assert gts[4] == (None, None)


def test_a_bad_mode_is_refused(tmp_path):
    with pytest.raises(SystemExit):
        _run(tmp_path, _SNP_WITH_LONE_STAR, singletons="stars")
