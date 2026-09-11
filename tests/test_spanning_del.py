"""`*` is a missingness annotation, not an allele.

A `*` in ALT says a deletion called somewhere else covers this position in some samples. It
names no base, and it does not say *which* deletion, so two samples carrying it cannot be
said to share an origin — which is the test that makes D384A worth separating from D384G.
Treating it as a third allele inflates every heterozygosity-based statistic with a
non-allele and makes an ordinary SNP look multiallelic.

The policy (MULTIALLELIC_PLAN.md §10, decision 4): drop the allele, keep the site, recode
the calls that name it as missing, and count what that cost.
"""

import shutil

import pytest

pysam = pytest.importorskip("pysam")
pytest.importorskip("cyvcf2")

from plasgenomicsutils.lib import vcf_filters as F
from plasgenomicsutils.lib.spanning_del import spanning_del_to_missing

needs = pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not on PATH")

_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n"
)


def _vcf(tmp_path, body, name="in.vcf"):
    p = tmp_path / name
    p.write_text(_HDR + body)
    return str(p)


def _rows(path):
    with pysam.VariantFile(path) as vf:
        return [(r.pos, r.alleles,
                 {s: r.samples[s].alleles for s in r.samples}) for r in vf]


def test_a_real_snp_under_a_deletion_survives_as_a_biallelic_snp(tmp_path):
    """`A > *,T` is the case the whole policy is for: 312 such records in the fixture."""
    src = _vcf(tmp_path, "chr1\t100\t.\tA\t*,T\t.\t.\t.\tGT:AD\t"
                         "0/0:9,0,0\t1/1:0,8,0\t2/2:0,0,7\t2/2:0,0,6\n")
    out = str(tmp_path / "o.vcf")
    st = F.spanning_del_filter(src, out)
    (pos, alleles, gts), = _rows(out)
    assert alleles == ("A", "T"), "the * allele is gone and T is re-indexed"
    assert gts["s1"] == ("A", "A")
    assert gts["s2"] == (None, None), "the * carrier is now missing"
    assert gts["s3"] == ("T", "T")
    assert st["records_with_spanning_del"] == 1
    assert st["calls_recoded"] == 1


def test_a_partial_call_keeps_the_allele_it_does_carry(tmp_path):
    """`*/T` means one haplotype is deleted and the other carries T.

    14.5% of the `*` calls in the shipped Pf7 fixture are partial like this. Discarding the
    whole call would throw away a real observation, so only the deleted slot goes missing.
    """
    src = _vcf(tmp_path, "chr1\t100\t.\tA\t*,T\t.\t.\t.\tGT:AD\t"
                         "0/0:9,0,0\t1/2:0,4,5\t2/2:0,0,7\t0/1:6,3,0\n")
    out = str(tmp_path / "o.vcf")
    st = F.spanning_del_filter(src, out)
    (_pos, alleles, gts), = _rows(out)
    assert alleles == ("A", "T")
    assert gts["s2"] == (None, "T"), "the T half of a */T call is kept"
    assert gts["s4"] == ("A", None), "and so is the A half of an A/* call"
    assert st["slots_recoded"] == 2
    assert st["calls_fully_missing"] == 0


def test_a_record_whose_only_alt_is_a_star_is_left_with_no_variant(tmp_path):
    """Half the polymorphic records in the fixture are this: not a SNP site at all."""
    src = _vcf(tmp_path, "chr1\t100\t.\tA\t*\t.\t.\t.\tGT:AD\t"
                         "0/0:9,0\t1/1:0,8\t1/1:0,7\t0/0:6,0\n")
    out = str(tmp_path / "o.vcf")
    st = F.spanning_del_filter(src, out)
    (_pos, alleles, _gts), = _rows(out)
    assert alleles == ("A",), "nothing but the reference is left"
    assert st["records_star_only"] == 1


def test_a_record_with_no_star_is_untouched(tmp_path):
    src = _vcf(tmp_path, "chr1\t100\t.\tA\tT,G\t.\t.\t.\tGT:AD\t"
                         "0/0:9,0,0\t1/1:0,8,0\t2/2:0,0,7\t1/2:0,4,4\n")
    out = str(tmp_path / "o.vcf")
    st = F.spanning_del_filter(src, out, trim=False)
    (_pos, alleles, gts), = _rows(out)
    assert alleles == ("A", "T", "G")
    assert gts["s4"] == ("T", "G")
    assert st["calls_recoded"] == 0
    assert st["records_with_spanning_del"] == 0


def test_the_depth_array_is_re_laid_with_the_alleles(tmp_path):
    """AD is Number=R, so dropping an allele has to drop its column too."""
    src = _vcf(tmp_path, "chr1\t100\t.\tA\t*,T\t.\t.\t.\tGT:AD\t"
                         "0/0:9,1,0\t1/1:0,8,0\t2/2:0,2,7\t2/2:0,0,6\n")
    out = str(tmp_path / "o.vcf")
    F.spanning_del_filter(src, out)
    with pysam.VariantFile(out) as vf:
        rec = next(iter(vf))
        assert len(rec.alleles) == 2
        assert tuple(rec.samples["s1"]["AD"]) == (9, 0), "the * column is gone, not shifted"
        assert tuple(rec.samples["s3"]["AD"]) == (0, 7)


def test_the_recode_alone_leaves_the_allele_in_place(tmp_path):
    """The two halves are separable: recoding is lossless about the ALT column."""
    src = _vcf(tmp_path, "chr1\t100\t.\tA\t*,T\t.\t.\t.\tGT:AD\t"
                         "0/0:9,0,0\t1/1:0,8,0\t2/2:0,0,7\t2/2:0,0,6\n")
    out = str(tmp_path / "o.vcf")
    st = spanning_del_to_missing(src, out)
    (_pos, alleles, gts), = _rows(out)
    assert alleles == ("A", "*", "T"), "recoding does not touch ALT; trimming does"
    assert gts["s2"] == (None, None)
    assert st["calls_recoded"] == 1


@needs
def test_the_policy_is_measured_on_the_shipped_fixture(tmp_path):
    """The scoping claim in the plan, checked against the data it was measured on."""
    import subprocess
    from pathlib import Path

    fixture = Path(__file__).parent / "data" / "ghana_cambodia.pf7.tiny.bcf"
    trimmed = str(tmp_path / "trim.bcf")
    subprocess.run(["bcftools", "view", "--trim-alt-alleles", str(fixture),
                    "-Ob", "-o", trimmed], check=True, stderr=subprocess.DEVNULL)

    def n_multiallelic(path):
        p = subprocess.run(["bcftools", "query", "-f", "%ALT\n", path],
                           stdout=subprocess.PIPE, text=True, stderr=subprocess.DEVNULL)
        return sum(1 for ln in p.stdout.split()
                   if len([a for a in ln.split(",") if a and a != "."]) >= 2)

    before = n_multiallelic(trimmed)
    out = str(tmp_path / "nostar.bcf")
    F.spanning_del_filter(trimmed, out)
    after = n_multiallelic(out)
    assert before == 544, f"the plan measured 544 multiallelic records, got {before}"
    assert after == 232, f"the plan measured 232 after the * policy, got {after}"


def test_the_step_is_in_the_default_chain_but_switched_off():
    """Off by default, and discoverable rather than absent.

    Recoding `*` to missing discards a real observation: a site with 20 deleted samples and
    5 carrying a variant then reads as though 20 samples could not be called there. In
    *P. falciparum* that is worse than usual, because a deletion across a dimorphic region
    is often the other haplotype rather than a dropout, and which samples carry it is a
    result rather than noise. So the recode is a choice the analysis makes, not a default
    the pipeline makes for it.
    """
    from plasgenomicsutils.lib.filter_pipeline import DEFAULT_CONFIG

    steps = {s["name"]: s for s in DEFAULT_CONFIG["steps"]}
    assert "spanning_del_filter" in steps, "it must stay visible in the written-out config"
    assert steps["spanning_del_filter"].get("enabled") is False

    # and where it sits still matters, because turning it on has to work
    names = [s["name"] for s in DEFAULT_CONFIG["steps"]]
    assert names.index("spanning_del_filter") < names.index("biallelic_snp_filter"), \
        "run after the biallelic test and the records it would rescue are already gone"


def test_turning_it_on_is_what_rescues_the_records(tmp_path):
    """The two orderings give different answers, which is why the position is pinned."""
    src = _vcf(tmp_path, "chr1\t100\t.\tA\t*,T\t.\t.\t.\tGT:AD\t"
                         "0/0:9,0,0\t1/1:0,8,0\t2/2:0,0,7\t2/2:0,0,6\n")
    # left alone, the record is multiallelic and a biallelic test removes it
    kept = str(tmp_path / "kept.bcf")
    F.biallelic_snp_filter(src, kept, trim=True, snps_only=True, biallelic=True,
                           mnp_handling="remove")
    with pysam.VariantFile(kept) as vf:
        assert len(list(vf)) == 0

    # recoded first, the same record survives as the biallelic SNP it always was
    recoded = str(tmp_path / "recoded.bcf")
    F.spanning_del_filter(src, recoded)
    out = str(tmp_path / "out.bcf")
    F.biallelic_snp_filter(recoded, out, trim=True, snps_only=True, biallelic=True,
                           mnp_handling="remove")
    with pysam.VariantFile(out) as vf:
        rec, = list(vf)
        assert rec.alleles == ("A", "T")
