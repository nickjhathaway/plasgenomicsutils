"""The per-allele singleton filter: recode the singleton alternate's calls and trim it,
instead of keeping the whole record because *some other* alternate is well supported.

The record-level drop (the default) is right when the question is "does this record hold a
supported variant". It is wrong when a record holds a good alternate AND a singleton one: the
record survives and the singleton allele rides along on it, into every count built on the
panel. `singleton_to_missing` blanks the calls that name such an alternate so a following
`--trim-alt-alleles` removes it. `*` is never touched here -- a spanning deletion is not a
variant allele, and a singleton one is `spanning_del_filter`'s business.
"""

import shutil

import pytest

pysam = pytest.importorskip("pysam")
pytest.importorskip("cyvcf2")

from plasgenomicsutils.lib.singletons import singleton_to_missing

_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
)

_SAMPLES = [f"s{i}" for i in range(1, 9)]


def _cells(alt, gts):
    n_ad = len(alt.split(",")) + 1
    out = []
    for gt in gts:
        k = int(gt.split("/")[0])
        ad = ["0"] * n_ad
        ad[k] = "40"
        out.append(f"{gt}:" + ",".join(ad))
    return out


def _vcf(tmp_path, rows, name="in.vcf"):
    hdr = _HDR + ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(_SAMPLES) + "\n")
    body = "".join(f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT:AD\t" + "\t".join(_cells(alt, gts)) + "\n"
                   for pos, alt, gts in rows)
    p = tmp_path / name
    p.write_text(hdr + body)
    return str(p)


def _recode(tmp_path, rows, min_samples=1):
    out = str(tmp_path / "recoded.vcf")
    st = singleton_to_missing(_vcf(tmp_path, rows), out, min_samples=min_samples)
    with pysam.VariantFile(out) as f:
        recs = [(r.pos, list(r.alts or ()),
                 [tuple(r.samples[s]["GT"]) for s in _SAMPLES]) for r in f]
    return recs, st


def test_a_singleton_alternate_beside_a_good_one_is_blanked_not_kept(tmp_path):
    # ALT1 (T) carried by 4 samples; ALT2 (G) private to one. The T calls stay; the G call
    # goes missing so a trim will remove G.
    rows = [(1000, "T,G", ["1/1", "1/1", "1/1", "1/1", "2/2", "0/0", "0/0", "0/0"])]
    recs, st = _recode(tmp_path, rows)
    (_pos, alts, gts), = recs
    assert alts == ["T", "G"]                       # ALT column untouched (trim is separate)
    assert gts[4] == (None, None)                   # the G carrier, blanked
    assert gts[0] == (1, 1)                          # a T carrier, untouched
    assert st["alts_recoded"] == 1
    assert st["records_with_singleton_alt"] == 1
    assert st["records_all_alts_singleton"] == 0


def test_a_record_whose_every_alternate_is_a_singleton_is_flagged(tmp_path):
    rows = [(1000, "T,G", ["1/1", "2/2", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"])]
    recs, st = _recode(tmp_path, rows)
    (_pos, _alts, gts), = recs
    assert gts[0] == (None, None) and gts[1] == (None, None)
    assert st["records_all_alts_singleton"] == 1


def test_a_star_allele_is_never_recoded(tmp_path):
    # `*` is not a variant allele; even a singleton one is left for spanning_del_filter
    rows = [(1000, "T,*", ["1/1", "1/1", "0/0", "2/2", "0/0", "0/0", "0/0", "0/0"])]
    recs, st = _recode(tmp_path, rows)
    (_pos, _alts, gts), = recs
    assert gts[3] == (2, 2)                           # the lone `*` carrier, untouched
    assert st["alts_recoded"] == 0
    assert st["records_with_singleton_alt"] == 0


def test_a_star_only_record_is_written_through(tmp_path):
    rows = [(1000, "*", ["1/1", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"])]
    recs, st = _recode(tmp_path, rows)
    (_pos, alts, gts), = recs
    assert alts == ["*"] and gts[0] == (1, 1)
    assert st["records_with_singleton_alt"] == 0


def test_min_samples_raises_the_bar(tmp_path):
    # with min_samples=2, an allele carried by exactly 2 is a singleton here too
    rows = [(1000, "T,G", ["1/1", "1/1", "2/2", "2/2", "2/2", "0/0", "0/0", "0/0"])]
    recs, st = _recode(tmp_path, rows, min_samples=2)
    (_pos, _alts, gts), = recs
    assert gts[0] == (None, None) and gts[1] == (None, None)   # T carried by 2: recoded
    assert gts[2] == (2, 2)                                     # G carried by 3: kept
    assert st["alts_recoded"] == 1


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_end_to_end_keeps_the_good_allele_and_drops_the_all_singleton_record(tmp_path):
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [
        (1000, "T,G", ["1/1", "1/1", "1/1", "1/1", "2/2", "0/0", "0/0", "0/0"]),  # T good, G singleton
        (2000, "T,G", ["1/1", "2/2", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"]),  # both singleton
        (3000, "C", ["1/1", "1/1", "1/1", "0/0", "0/0", "0/0", "0/0", "0/0"]),    # good biallelic
    ]
    src = _vcf(tmp_path, rows)
    out = str(tmp_path / "filtered.vcf")
    singleton_add_ads(src, out, min_samples=1, per_allele=True)
    with pysam.VariantFile(out) as f:
        kept = [(r.pos, list(r.alts or ())) for r in f]
    # 1000 survives with only T (G trimmed); 2000 gone; 3000 survives
    assert (1000, ["T"]) in kept
    assert (3000, ["C"]) in kept
    assert all(pos != 2000 for pos, _ in kept)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_default_still_keeps_the_whole_record(tmp_path):
    # per_allele defaults off: the record with a good allele survives, singleton allele and all
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [(1000, "T,G", ["1/1", "1/1", "1/1", "1/1", "2/2", "0/0", "0/0", "0/0"])]
    src = _vcf(tmp_path, rows)
    out = str(tmp_path / "filtered.vcf")
    singleton_add_ads(src, out, min_samples=1)
    with pysam.VariantFile(out) as f:
        kept = [(r.pos, list(r.alts or ())) for r in f]
    assert kept == [(1000, ["T", "G"])]


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_a_singleton_snp_beside_a_common_deletion_leaves_the_deletion(tmp_path):
    """The subtle case, and the reason per-allele mode can keep MORE records than the
    record-level default. `REF=A ALT=T,*` with T private to one sample and `*` common: the
    singleton SNP is blanked and trimmed away, and what remains is a real spanning deletion.
    Per-allele mode keeps it (as a `*`-only record, for spanning_del_filter to judge); the
    record-level default drops the whole site because its only real allele is a singleton,
    throwing the deletion away with it."""
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [(1000, "T,*", ["1/1"] + ["2/2"] * 6 + ["0/0"])]   # T: 1 carrier, *: 6 carriers
    src = _vcf(tmp_path, rows)

    per = str(tmp_path / "per.vcf")
    singleton_add_ads(src, per, min_samples=1, per_allele=True)
    with pysam.VariantFile(per) as f:
        recs = [(r.pos, list(r.alts or ())) for r in f]
    assert recs == [(1000, ["*"])]              # the singleton SNP gone, the deletion kept

    rec = str(tmp_path / "rec.vcf")
    singleton_add_ads(src, rec, min_samples=1)  # default: whole site dropped
    with pysam.VariantFile(rec) as f:
        assert [r.pos for r in f] == []


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_a_singleton_snp_with_no_deletion_is_dropped_outright(tmp_path):
    """No `*` to fall back on: the record is ref-only after the trim and dropped."""
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [(1000, "T", ["1/1", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"])]
    src = _vcf(tmp_path, rows)
    out = str(tmp_path / "out.vcf")
    singleton_add_ads(src, out, min_samples=1, per_allele=True)
    with pysam.VariantFile(out) as f:
        assert [r.pos for r in f] == []
