"""`*` is not a variant allele, and a singleton filter that counts it drops the wrong records.

A spanning deletion says the sequence is not there on that haplotype. That is a real
observation -- it is why `spanning_del_filter` is default-off -- but it is not one of the
alternates a SNP filter is about. `MAX(INFO/AC_SAMP)` cannot tell the two apart, so on a
record like `REF=A ALT=T,*` where T is private to one sample and `*` is carried by thirty,
the maximum is 30 and the record survives a filter whose entire purpose is to drop records
with no well-supported alternate. The SNP it claims to hold is a singleton.

`AC_SAMP_MAX` is the scalar that answers the question actually being asked: the best carrier
count over the *real* alternates, falling back to the star's own count when `*` is the only
alternate -- there the deletion IS the variant, and dropping a well-attested one as a
singleton would be a different rule than the one being applied.
"""

import shutil

import pytest

pysam = pytest.importorskip("pysam")

from plasgenomicsutils.lib.allele_counts import (
    ALT_SAMPLE_MAX_TAG, ALT_SAMPLE_TAG, add_alt_sample_counts,
)

_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
)

_SAMPLES = [f"s{i}" for i in range(1, 9)]


def _vcf(tmp_path, rows, name="in.vcf"):
    hdr = _HDR + ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(_SAMPLES) + "\n")
    def _cells(alt, gts):
        n_ad = len(alt.split(",")) + 1
        out = []
        for gt in gts:
            k = int(gt.split("/")[0])
            ad = ["0"] * n_ad
            ad[k] = "40"
            out.append(f"{gt}:" + ",".join(ad))
        return out

    body = "".join(f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT:AD\t"
                   + "\t".join(_cells(alt, gts)) + "\n"
                   for pos, alt, gts in rows)
    p = tmp_path / name
    p.write_text(hdr + body)
    return str(p)


def _tagged(tmp_path, rows):
    out = str(tmp_path / "out.vcf")
    st = add_alt_sample_counts(_vcf(tmp_path, rows), out)
    with pysam.VariantFile(out) as f:
        recs = [(r.pos, tuple(r.info[ALT_SAMPLE_TAG]), r.info[ALT_SAMPLE_MAX_TAG])
                for r in f]
    return recs, st


def test_a_private_snp_beside_a_common_star_is_still_a_singleton(tmp_path):
    # one carrier of T, six of `*`: MAX(AC_SAMP) says 6, and the record holds no real allele
    # carried by more than one sample
    rows = [(1000, "T,*", ["1/1"] + ["2/2"] * 6 + ["0/0"])]
    recs, st = _tagged(tmp_path, rows)
    (_pos, ac, best), = recs
    assert ac == (1, 6)          # AC_SAMP stays faithful: `*` is counted, as Number=A needs
    assert max(ac) == 6
    assert best == 1             # ...and the scalar the filter reads does not
    assert st["records_star_only_support"] == 1


def test_a_star_only_record_keeps_its_own_count(tmp_path):
    # nothing to fall back to, and here the deletion IS the variant -- dropping a deletion
    # carried by six samples as a "singleton" would be a different rule
    rows = [(2000, "*", ["1/1"] * 6 + ["0/0"] * 2)]
    recs, _st = _tagged(tmp_path, rows)
    (_pos, ac, best), = recs
    assert ac == (6,)
    assert best == 6


def test_a_record_with_no_star_is_unchanged(tmp_path):
    rows = [(3000, "T,G", ["1/1"] * 3 + ["2/2"] * 4 + ["0/0"])]
    recs, st = _tagged(tmp_path, rows)
    (_pos, ac, best), = recs
    assert ac == (3, 4)
    assert best == max(ac) == 4
    assert st["records_with_star"] == 0


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_the_singleton_filter_drops_it_end_to_end(tmp_path):
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [
        (1000, "T,*", ["1/1"] + ["2/2"] * 6 + ["0/0"]),   # no real allele above 1 carrier
        (2000, "T,*", ["1/1", "1/1"] + ["2/2"] * 5 + ["0/0"]),  # T has 2: kept
        (3000, "*", ["1/1"] * 6 + ["0/0"] * 2),          # deletion is the variant: kept
    ]
    src = _vcf(tmp_path, rows)
    out = str(tmp_path / "filtered.vcf")
    singleton_add_ads(src, out, min_samples=1)
    with pysam.VariantFile(out) as f:
        kept = sorted(r.pos for r in f)
    assert kept == [2000, 3000]
