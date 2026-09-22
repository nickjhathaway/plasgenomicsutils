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


# --- singleton `*` ----------------------------------------------------------------------
#
# A `*` is never *support* for a record, but a `*` carried by one sample is a private
# observation like any other, and left on the record it takes the site out of any
# `--snps-only` panel. So the singleton filter blanks it, in both modes; a `*` with more
# carriers is left for `spanning_del_filter`.


def test_the_default_mode_leaves_real_alleles_alone_but_blanks_a_singleton_star(tmp_path):
    rows = [(1000, "T,*", ["1/1"] * 3 + ["2/2"] + ["0/0"] * 4)]     # T: 3 carriers, *: 1
    out = str(tmp_path / "star.vcf")
    st = singleton_to_missing(_vcf(tmp_path, rows), out, min_samples=1, alleles="star")
    with pysam.VariantFile(out) as f:
        (rec,) = list(f)
        gts = [tuple(rec.samples[s]["GT"]) for s in _SAMPLES]
    assert list(rec.alts) == ["T", "*"]                     # ALT untouched: the trim's job
    assert gts[3] == (None, None)                           # the one `*` carrier is blanked
    assert gts[:3] == [(1, 1)] * 3 and gts[4:] == [(0, 0)] * 4
    assert st["stars_recoded"] == 1 and st["records_with_singleton_star"] == 1
    assert st["alts_recoded"] == 0 and st["records_with_singleton_alt"] == 0


def test_star_mode_never_touches_a_singleton_real_allele(tmp_path):
    rows = [(1000, "T,*", ["1/1"] + ["2/2"] * 3 + ["0/0"] * 4)]     # T: 1 carrier, *: 3
    out = str(tmp_path / "star.vcf")
    st = singleton_to_missing(_vcf(tmp_path, rows), out, min_samples=1, alleles="star")
    with pysam.VariantFile(out) as f:
        (rec,) = list(f)
        gts = [tuple(rec.samples[s]["GT"]) for s in _SAMPLES]
    assert gts[0] == (1, 1)                                 # the singleton SNP is not its job
    assert st["stars_recoded"] == 0 and st["alts_recoded"] == 0


def test_all_mode_blanks_both(tmp_path):
    rows = [(1000, "T,G,*", ["1/1"] * 3 + ["2/2"] + ["3/3"] + ["0/0"] * 3)]
    out = str(tmp_path / "all.vcf")
    st = singleton_to_missing(_vcf(tmp_path, rows), out, min_samples=1, alleles="all")
    with pysam.VariantFile(out) as f:
        (rec,) = list(f)
        gts = [tuple(rec.samples[s]["GT"]) for s in _SAMPLES]
    assert gts[3] == (None, None) and gts[4] == (None, None)
    assert gts[:3] == [(1, 1)] * 3
    assert st["alts_recoded"] == 1 and st["stars_recoded"] == 1
    assert st["records_all_alts_singleton"] == 0            # T survives, so not all-singleton


def test_a_bad_mode_is_refused(tmp_path):
    with pytest.raises(ValueError):
        singleton_to_missing(_vcf(tmp_path, []), str(tmp_path / "x.vcf"), alleles="stars")


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_a_common_snp_under_one_samples_deletion_comes_out_as_a_plain_snp(tmp_path):
    """The case this exists for: `A > T,*` with T in three samples and `*` in one. Before,
    the record-level filter kept the record with its `*` and `--snps-only` then threw the
    whole site away. Now the one deleted sample is blanked, the `*` trimmed, and the record
    is the biallelic SNP it always was for everyone else -- in both modes."""
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [(1000, "T,*", ["1/1"] * 3 + ["2/2"] + ["0/0"] * 4)]
    src = _vcf(tmp_path, rows)
    for per_allele in (False, True):
        out = str(tmp_path / f"out_{per_allele}.vcf")
        singleton_add_ads(src, out, min_samples=1, per_allele=per_allele)
        with pysam.VariantFile(out) as f:
            (rec,) = list(f)
            gts = [tuple(rec.samples[s]["GT"]) for s in _SAMPLES]
        assert list(rec.alts) == ["T"], per_allele
        assert gts[3] == (None, None), per_allele
        assert gts[:3] == [(1, 1)] * 3, per_allele


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_a_common_deletion_is_left_for_spanning_del_filter(tmp_path):
    """Two `*` carriers with min_samples=1 is not a singleton: the deletion stays, and the
    record keeps its `*` for `spanning_del_filter` (or `--snps-only`) to judge."""
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [(1000, "T,*", ["1/1"] * 3 + ["2/2"] * 2 + ["0/0"] * 3)]
    src = _vcf(tmp_path, rows)
    out = str(tmp_path / "out.vcf")
    singleton_add_ads(src, out, min_samples=1)
    with pysam.VariantFile(out) as f:
        (rec,) = list(f)
        gts = [tuple(rec.samples[s]["GT"]) for s in _SAMPLES]
    assert list(rec.alts) == ["T", "*"]
    assert gts[3] == (2, 2) and gts[4] == (2, 2)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_a_singleton_star_only_record_is_dropped_and_a_common_one_kept(tmp_path):
    """A `*`-only record is judged on the star's own count. A singleton one is blanked to
    ref-only and dropped; a common one is untouched and kept as before."""
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [(1000, "*", ["1/1"] + ["0/0"] * 7),
            (2000, "*", ["1/1"] * 4 + ["0/0"] * 4)]
    src = _vcf(tmp_path, rows)
    out = str(tmp_path / "out.vcf")
    singleton_add_ads(src, out, min_samples=1)
    with pysam.VariantFile(out) as f:
        kept = [(r.pos, list(r.alts or ())) for r in f]
    assert kept == [(2000, ["*"])]


def test_a_star_only_record_is_never_recoded_however_private(tmp_path):
    """No real alternate means no SNP being held hostage: the deletion IS the variant here,
    the same carve-out `AC_SAMP_MAX` makes. Recoding would empty the ALT column, and a
    whitelist rescue would then put an allele-less record back instead of the deletion."""
    rows = [(1000, "*", ["1/1"] + ["0/0"] * 7)]
    out = str(tmp_path / "star_only.vcf")
    st = singleton_to_missing(_vcf(tmp_path, rows), out, min_samples=1, alleles="all")
    with pysam.VariantFile(out) as f:
        (rec,) = list(f)
        gts = [tuple(rec.samples[s]["GT"]) for s in _SAMPLES]
    assert list(rec.alts) == ["*"] and gts[0] == (1, 1)
    assert st["stars_recoded"] == 0


def test_a_star_beside_an_uncarried_real_allele_is_left_alone(tmp_path):
    """The general form of the same rule. `T` is in the ALT column but nobody carries it --
    a leftover from an earlier step's re-genotyping -- so blanking the `*` would rescue no
    SNP and would leave the record with no allele at all once trimmed."""
    rows = [(1000, "T,*", ["2/2"] + ["0/0"] * 7)]      # T: 0 carriers, *: 1
    out = str(tmp_path / "uncarried.vcf")
    st = singleton_to_missing(_vcf(tmp_path, rows), out, min_samples=1, alleles="star")
    with pysam.VariantFile(out) as f:
        (rec,) = list(f)
        gts = [tuple(rec.samples[s]["GT"]) for s in _SAMPLES]
    assert list(rec.alts) == ["T", "*"] and gts[0] == (2, 2)
    assert st["stars_recoded"] == 0


def test_a_star_is_still_recoded_when_the_other_singleton_alt_is_not_its_only_company(tmp_path):
    """`all` mode, three alternates: G is a singleton and goes, T is well supported and
    stays, so the `*` has a SNP to rescue and goes too."""
    rows = [(1000, "T,G,*", ["1/1"] * 3 + ["2/2"] + ["3/3"] + ["0/0"] * 3)]
    out = str(tmp_path / "three.vcf")
    st = singleton_to_missing(_vcf(tmp_path, rows), out, min_samples=1, alleles="all")
    assert st["stars_recoded"] == 1 and st["alts_recoded"] == 1


def test_a_star_is_not_recoded_when_every_real_alt_is_also_a_singleton(tmp_path):
    """`all` mode again, but nothing real survives: T is a singleton and goes, leaving the
    `*` as the record's only content. It stays, so the record is the deletion it is rather
    than an empty one -- which is what `spanning_del_filter` is there to judge."""
    rows = [(1000, "T,*", ["1/1"] + ["2/2"] + ["0/0"] * 6)]
    out = str(tmp_path / "both.vcf")
    st = singleton_to_missing(_vcf(tmp_path, rows), out, min_samples=1, alleles="all")
    with pysam.VariantFile(out) as f:
        (rec,) = list(f)
        gts = [tuple(rec.samples[s]["GT"]) for s in _SAMPLES]
    assert gts[0] == (None, None) and gts[1] == (2, 2)
    assert st["alts_recoded"] == 1 and st["stars_recoded"] == 0


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")
def test_a_whitelisted_record_is_exempt_from_the_drop_not_from_the_recode(tmp_path):
    """What the whitelist promises and what it does not. The singleton record comes back
    whole (exempt from this filter's drop). The record beside a singleton `*` is still
    recoded -- the whitelist is not a freeze on a record's genotypes -- but what it keeps is
    a clean SNP, which is what naming the position asked for. A `*`-only record is rescued as
    the deletion it is, not as an empty ALT."""
    from plasgenomicsutils.lib.vcf_filters import singleton_add_ads

    rows = [(1000, "T,*", ["1/1"] * 3 + ["2/2"] + ["0/0"] * 4),   # common SNP, singleton `*`
            (2000, "T", ["1/1"] + ["0/0"] * 7),                   # singleton SNP: dropped
            (3000, "*", ["1/1"] + ["0/0"] * 7)]                   # singleton `*`-only
    src = _vcf(tmp_path, rows)
    bed = tmp_path / "wl.bed"
    bed.write_text("chr1\t999\t1000\nchr1\t1999\t2000\nchr1\t2999\t3000\n")
    out = str(tmp_path / "out.vcf")
    singleton_add_ads(src, out, min_samples=1, keep_bed=str(bed))
    with pysam.VariantFile(out) as f:
        got = {r.pos: (list(r.alts or []),
                       [tuple(r.samples[s]["GT"]) for s in _SAMPLES]) for r in f}
    assert got[1000][0] == ["T"]                       # `*` trimmed off
    assert got[1000][1][3] == (None, None)             # its one carrier blanked
    assert got[1000][1][:3] == [(1, 1)] * 3            # everyone else untouched
    assert got[2000][0] == ["T"] and got[2000][1][0] == (1, 1)   # rescued whole
    assert got[3000][0] == ["*"] and got[3000][1][0] == (1, 1)   # rescued as the deletion
