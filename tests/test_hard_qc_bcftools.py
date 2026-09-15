"""hard_qc_filter --caller bcftools: the same questions asked of a bcftools callset.

GATK's metrics are not in a bcftools VCF, so the filter reads what `bcftools mpileup`
writes instead. Two of the mappings are not straight renames, and both are covered here:
bcftools FS is the strand-bias p-value rather than a Phred score, and the *BZ tags are
two-sided.
"""

import math
import pathlib
import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import vcf_filters as F

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"),
                                reason="bcftools not on PATH")

BCF_INFO = {"DP": ("1", "Integer"), "MQ": ("1", "Integer"), "FS": ("1", "Float"),
            "RPBZ": ("1", "Float"), "SCBZ": ("1", "Float"), "MQBZ": ("1", "Float"),
            "MQSBZ": ("1", "Float"), "BQBZ": ("1", "Float"), "MQ0F": ("1", "Float"),
            "ADF": ("R", "Integer"), "ADR": ("R", "Integer")}


def _vcf(tmp_path, records, tags=tuple(BCF_INFO), name="in.vcf"):
    """A minimal bcftools-shaped callset; `records` are (pos, {INFO tag: value})."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>"]
    for t in tags:
        n, ty = BCF_INFO[t]
        hdr.append(f'##INFO=<ID={t},Number={n},Type={ty},Description="{t}">')
    hdr.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
    hdr.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1")
    for pos, info in records:
        kv = ";".join(f"{k}={v}" for k, v in info.items())
        hdr.append(f"chr1\t{pos}\t.\tA\tT\t222\t.\t{kv}\tGT\t0/1")
    p = tmp_path / name
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _kept(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\n", path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return [int(x) for x in out.stdout.split()]


# a clean site: balanced strands on both alleles, every bias statistic at rest
CLEAN = {"DP": 40, "MQ": 60, "FS": 0.25, "RPBZ": 0.0, "SCBZ": 0.0, "MQBZ": 0.0,
         "MQSBZ": 0.0, "BQBZ": 0.0, "MQ0F": 0.0, "ADF": "10,10", "ADR": "10,10"}


def _strands(ref_f, ref_r, alt_f, alt_r):
    """INFO overrides putting a given 2x2 strand table on a record."""
    return {"ADF": f"{ref_f},{alt_f}", "ADR": f"{ref_r},{alt_r}"}


def _sor(ref_f, ref_r, alt_f, alt_r):
    """GATK's StrandOddsRatio, written out here so the test does not reuse the filter's
    own arithmetic to check the filter's own arithmetic."""
    a, b, c, d = ref_f + 1, ref_r + 1, alt_f + 1, alt_r + 1
    ratio = (a / b) * (d / c)
    return (math.log(ratio + 1 / ratio) + math.log(min(a, b) / max(a, b))
            - math.log(min(c, d) / max(c, d)))


def _with(**kw):
    d = dict(CLEAN)
    d.update(kw)
    return d


# ---- the two metrics the filter exists for ---------------------------------------

def test_strand_bias_is_a_p_value_so_small_is_bad(tmp_path):
    # bcftools FS is Fisher's p-value, not GATK's Phred score: a biased site has a tiny
    # FS. Filtering it the GATK way round (FS > 60) would keep exactly the wrong records.
    # It is opt-in -- see the strand-bias section below for why -- but when asked for, it
    # has to be asked the right way round.
    vcf = _vcf(tmp_path, [(1000, CLEAN), (2000, _with(FS=7.2e-12))])
    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools", strand_bias_p=1e-6)
    assert _kept(out) == [1000]


def test_variants_at_the_ends_of_reads_are_dropped_either_way_round(tmp_path):
    # RPBZ is two-sided ("closer to 0 is better"), unlike GATK's one-sided
    # ReadPosRankSum < -5, so a positive z has to be caught as well as a negative one.
    vcf = _vcf(tmp_path, [(1000, CLEAN), (2000, _with(RPBZ=6.2)), (3000, _with(RPBZ=-6.2))])
    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools")
    assert _kept(out) == [1000]


def test_soft_clipping_counts_as_a_read_end_signal(tmp_path):
    # SCBZ rides the same threshold: reads clipped at their ends are the other way the
    # same artifact shows up
    vcf = _vcf(tmp_path, [(1000, CLEAN), (2000, _with(SCBZ=-8.0))])
    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools")
    assert _kept(out) == [1000]


def test_each_threshold_can_be_switched_off(tmp_path):
    biased = [(1000, CLEAN), (2000, _with(FS=1e-12)), (3000, _with(RPBZ=9.0))]
    vcf = _vcf(tmp_path, biased)
    out = str(tmp_path / "a.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools", strand_bias_p=None)
    assert _kept(out) == [1000, 2000]
    out = str(tmp_path / "b.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools", read_pos_z=None, strand_bias_p=1e-6)
    assert _kept(out) == [1000, 3000]


def test_mapping_quality_bias_and_the_optional_extras(tmp_path):
    vcf = _vcf(tmp_path, [(1000, CLEAN), (2000, _with(MQBZ=-7.0)), (3000, _with(MQSBZ=7.0)),
                          (4000, _with(BQBZ=9.0)), (5000, _with(MQ0F=0.5))])
    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools")           # BQBZ / MQ0F off by default
    assert _kept(out) == [1000, 4000, 5000]
    out2 = str(tmp_path / "out2.vcf")
    F.hard_qc_filter(vcf, out2, caller="bcftools", bqbz_z=5.0, mq0f=0.1)
    assert _kept(out2) == [1000]


# ---- QD does not carry across ----------------------------------------------------

def test_qd_is_off_by_default_for_bcftools_but_10_for_gatk(tmp_path):
    # QUAL 222 over DP 40 is 5.6: GATK-mode's QD >= 10 would discard a good bcftools callset,
    # so "auto" leaves it off there
    vcf = _vcf(tmp_path, [(1000, CLEAN)])
    out = str(tmp_path / "auto.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools")
    assert _kept(out) == [1000]
    # asking for it explicitly still works, and does drop the record
    out2 = str(tmp_path / "qd.vcf")
    F.hard_qc_filter(vcf, out2, caller="bcftools", qd=10)
    assert _kept(out2) == []
    # and a modest threshold on the right scale keeps it
    out3 = str(tmp_path / "qd2.vcf")
    F.hard_qc_filter(vcf, out3, caller="bcftools", qd=2)
    assert _kept(out3) == [1000]


def test_depth_resolves_to_info_dp(tmp_path):
    """A bare `DP` in a bcftools expression means FORMAT/DP where both exist."""
    assert "QUAL/INFO/DP" in F._bcftools_qc_expr(
        qd=20, mq=None, sor=None, strand_bias_p=None, max_bias_z=None, read_pos_z=None,
        mq0f=None, bqbz_z=None, bias_eff=None)


# ---- refusing to run rather than silently keeping everything ----------------------

def test_a_callset_without_the_tags_is_an_error_naming_them(tmp_path):
    # plain `bcftools mpileup | call` emits RPBZ and MQBZ but not FS, and a comparison
    # against an absent tag is false -- so the filter would quietly keep everything
    vcf = _vcf(tmp_path, [(1000, {"DP": 40, "MQ": 60, "RPBZ": 0.0})],
               tags=("DP", "MQ", "RPBZ"))
    with pytest.raises(SystemExit) as e:
        F.hard_qc_filter(vcf, str(tmp_path / "out.vcf"), caller="bcftools")
    msg = str(e.value)
    assert "INFO/SCBZ" in msg and "INFO/ADF" in msg
    assert "bcftools mpileup" in msg          # the call that would produce them
    # ...and disabling the thresholds that read the missing tags is enough to proceed
    out = str(tmp_path / "ok.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools", sor=None, read_pos_z=None,
                     max_bias_z=None)
    assert _kept(out) == [1000]


def test_an_unknown_caller_is_rejected(tmp_path):
    vcf = _vcf(tmp_path, [(1000, CLEAN)])
    with pytest.raises(SystemExit, match="caller must be gatk or bcftools"):
        F.hard_qc_filter(vcf, str(tmp_path / "o.vcf"), caller="freebayes")


def test_disabling_everything_is_an_error_not_a_no_op(tmp_path):
    vcf = _vcf(tmp_path, [(1000, CLEAN)])
    with pytest.raises(SystemExit, match="every bcftools threshold was disabled"):
        F.hard_qc_filter(vcf, str(tmp_path / "o.vcf"), caller="bcftools", qd=None,
                         mq=None, sor=None, strand_bias_p=None, read_pos_z=None,
                         max_bias_z=None)


# ---- the GATK path is untouched --------------------------------------------------

def test_gatk_mode_still_reads_gatk_metrics(tmp_path):
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>"]
    for t in ("QD", "MQ", "SOR", "MQRankSum", "ReadPosRankSum"):
        hdr.append(f'##INFO=<ID={t},Number=1,Type=Float,Description="{t}">')
    hdr.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
    hdr.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1")
    good = "QD=30;MQ=60;SOR=1;MQRankSum=0;ReadPosRankSum=0"
    bad = "QD=5;MQ=60;SOR=1;MQRankSum=0;ReadPosRankSum=0"
    hdr.append(f"chr1\t1000\t.\tA\tT\t222\t.\t{good}\tGT\t0/1")
    hdr.append(f"chr1\t2000\t.\tA\tT\t222\t.\t{bad}\tGT\t0/1")
    vcf = tmp_path / "gatk.vcf"
    vcf.write_text("\n".join(hdr) + "\n")

    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(str(vcf), out)                    # caller defaults to gatk
    assert _kept(out) == [1000]


# ---- non-variant records, counted on their own ------------------------------------
# The reason this is its own step: bcftools computes the bias statistics whether or not
# an ALT was called, so a hard QC rule removes non-variant records too, and the count of
# "nothing to call here" gets mixed into "failed quality".

def _mixed_vcf(tmp_path):
    """Two variant records and two non-variant ones, one of each biased."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>"]
    for t, (n, ty) in BCF_INFO.items():
        hdr.append(f'##INFO=<ID={t},Number={n},Type={ty},Description="{t}">')
    hdr.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
    hdr.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1")
    solo = {"ADF": "10", "ADR": "10"}        # Number=R: one entry where there is no ALT
    rows = [(1000, "T", CLEAN), (2000, "T", _with(**_strands(100, 100, 20, 0))),
            (3000, ".", _with(**solo)), (4000, ".", _with(MQ=30, **solo))]
    for pos, alt, info in rows:
        kv = ";".join(f"{k}={v}" for k, v in info.items())
        hdr.append(f"chr1\t{pos}\t.\tA\t{alt}\t222\t.\t{kv}\tGT\t0/1")
    p = tmp_path / "mixed.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def test_a_hard_qc_rule_does_remove_non_variant_records(tmp_path):
    """Why the separate step exists: this is what happens without it."""
    out = str(tmp_path / "qc.vcf")
    F.hard_qc_filter(_mixed_vcf(tmp_path), out, caller="bcftools")
    # 4000 is non-variant AND badly mapped, so QC takes it out along with 2000. (The z
    # tests are the exception since v0.2.4: a no-ALT record has no alt reads to size the
    # effect on, so they leave it alone -- see test_a_non_variant_record_is_not_judged_on_z.)
    assert _kept(out) == [1000, 3000]


def test_no_alt_filter_removes_them_and_counts_them(tmp_path, capsys):
    out = str(tmp_path / "out.vcf")
    F.no_alt_filter(_mixed_vcf(tmp_path), out)
    assert _kept(out) == [1000, 2000]
    assert "2 record(s) have no ALT allele (dropped)" in capsys.readouterr().err


def test_keep_passes_them_through_but_still_counts(tmp_path, capsys):
    out = str(tmp_path / "out.vcf")
    F.no_alt_filter(_mixed_vcf(tmp_path), out, keep=True)
    assert _kept(out) == [1000, 2000, 3000, 4000]
    assert "2 record(s) have no ALT allele (kept)" in capsys.readouterr().err


def test_the_two_reasons_are_separable(tmp_path):
    """Running it first splits 'nothing to call' from 'failed quality'."""
    vcf = _mixed_vcf(tmp_path)
    stage1 = str(tmp_path / "1.vcf")
    F.no_alt_filter(vcf, stage1)
    stage2 = str(tmp_path / "2.vcf")
    F.hard_qc_filter(stage1, stage2, caller="bcftools")
    assert _kept(stage1) == [1000, 2000]     # 2 non-variant removed here
    assert _kept(stage2) == [1000]           # 1 real variant failed QC there


# ---- strand bias as an effect size, not a p-value ---------------------------------
#
# FS is a Fisher p-value pooled over every sample, so a fixed cutoff asks "am I sure there
# is a skew" -- and the answer turns yes for any skew at all once a cohort is large enough.
# SOR asks how big the skew is, which is what the filter is supposed to mean.

@pytest.mark.parametrize("table,biased", [
    ((100, 100, 20, 20), False),   # balanced ref, balanced alt
    ((100, 100, 18, 2), True),     # alt 90% on one strand
    ((100, 100, 20, 0), True),     # alt entirely on one strand
    ((100, 1, 20, 1), False),      # ref AND alt skewed alike: a lopsided pileup, not bias
    ((100, 1, 10, 10), False),     # alt balanced where the ref is not
])
def test_sor_drops_a_one_strand_alt_and_keeps_a_lopsided_pileup(tmp_path, table, biased):
    v = _vcf(tmp_path, [(1000, _with(**_strands(*table)))])
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(v, out, caller="bcftools")
    assert (_kept(out) == []) is biased
    assert (_sor(*table) > 3) is biased          # and the reference agrees on why


@pytest.mark.parametrize("scale", [1, 5, 20, 100, 5000])
def test_the_strand_bias_verdict_does_not_move_with_cohort_size(tmp_path, scale):
    """The whole point. Multiply every count by the same factor -- the same strand
    composition seen in a bigger cohort -- and the verdict has to be unchanged. A p-value
    threshold fails this test: its verdict flips to 'biased' as the counts grow."""
    ok, bad = (60, 40, 12, 8), (60, 40, 19, 1)
    v = _vcf(tmp_path, [(1000, _with(**_strands(*[x * scale for x in ok]))),
                        (2000, _with(**_strands(*[x * scale for x in bad])))])
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(v, out, caller="bcftools")
    assert _kept(out) == [1000]


def test_fs_is_off_by_default_and_can_be_switched_back_on(tmp_path):
    """A vanishing FS on a site with no real skew is what a large cohort produces; the
    default must not drop it, and asking for FS explicitly must still work."""
    rec = [(1000, _with(FS=1e-30, **_strands(60, 40, 12, 8)))]
    a, b = str(tmp_path / "a.bcf"), str(tmp_path / "b.bcf")
    F.hard_qc_filter(_vcf(tmp_path, rec), a, caller="bcftools")
    assert _kept(a) == [1000]
    F.hard_qc_filter(_vcf(tmp_path, rec, name="in2.vcf"), b, caller="bcftools",
                     strand_bias_p=1e-6)
    assert _kept(b) == []


def test_sor_can_be_switched_off(tmp_path):
    v = _vcf(tmp_path, [(1000, _with(**_strands(100, 100, 20, 0)))])
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(v, out, caller="bcftools", sor=None, strand_bias_p=1e-6)
    assert _kept(out) == [1000]


def test_sor_needs_the_strand_tags_and_the_error_names_the_flag(tmp_path):
    v = _vcf(tmp_path, [(1000, {k: v for k, v in CLEAN.items()
                                if k not in ("ADF", "ADR")})],
             tags=[t for t in BCF_INFO if t not in ("ADF", "ADR")])
    with pytest.raises(SystemExit) as e:
        F.hard_qc_filter(v, str(tmp_path / "o.bcf"), caller="bcftools")
    assert "INFO/ADF" in str(e.value) and "INFO/ADR" in str(e.value)
    assert "--sor" in str(e.value)


def test_a_non_variant_record_is_not_judged_on_strand_bias(tmp_path):
    """ADF/ADR have no ALT entry on a no-ALT record, and an out-of-range index is not a
    verdict. no_alt_filter is what removes those, in its own step."""
    hdr = _vcf(tmp_path, [(1000, _with())])
    txt = pathlib.Path(hdr).read_text().replace(
        "chr1\t1000\t.\tA\tT\t222\t.\tDP=40",
        "chr1\t1000\t.\tA\t.\t222\t.\tDP=40").replace("ADF=10,10", "ADF=10").replace(
        "ADR=10,10", "ADR=10")
    v = tmp_path / "noalt.vcf"
    v.write_text(txt)
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(str(v), out, caller="bcftools")
    assert _kept(out) == [1000]


# ---- the z tests as significant *and* large ------------------------------------------
#
# RPBZ/SCBZ/MQBZ/MQSBZ are Mann-Whitney z-scores pooled over every read at the site, so a
# fixed z cutoff tightens as a cohort grows exactly as FS does. The effect size behind the
# z does not, so a z only counts when its effect is at least `bias_eff`. The counts come
# from ADF/ADR: ref vs alt for the allele tags, forward vs reverse for MQSBZ.

def _eff(z, n1, n2):
    """The effect size the filter is meant to recover, written out independently."""
    return abs(z) * math.sqrt((n1 + n2 + 1) / (12 * n1 * n2))


def _z_for(eff, n1, n2):
    """The z a given effect scores with these many reads behind it."""
    return eff / math.sqrt((n1 + n2 + 1) / (12 * n1 * n2))


@pytest.mark.parametrize("scale", [1, 10, 100, 1000])
def test_the_read_position_verdict_does_not_move_with_cohort_size(tmp_path, scale):
    """The whole point, as for SOR. The same modest shift (eff 0.05) and the same real
    artifact (eff 0.4), each seen in a cohort `scale` times larger: the z grows with
    sqrt(scale) but the verdict must not change."""
    n1, n2 = 100 * scale, 100 * scale
    counts = _strands(n1 // 2, n1 - n1 // 2, n2 // 2, n2 - n2 // 2)
    modest, artifact = _z_for(0.05, n1, n2), _z_for(0.4, n1, n2)
    v = _vcf(tmp_path, [(1000, _with(RPBZ=modest, **counts)),
                        (2000, _with(RPBZ=-artifact, **counts))])
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(v, out, caller="bcftools")
    assert _kept(out) == [1000]
    # ...and the plain z test, which is what bias_eff=None restores, does flip on the
    # modest shift once the cohort is large: z = 24 at scale 1000
    out2 = str(tmp_path / "o2.bcf")
    F.hard_qc_filter(v, out2, caller="bcftools", bias_eff=None)
    assert _kept(out2) == ([1000] if modest <= 5 else [])


def test_a_real_cohort_site_z_32_with_a_65_35_split_is_kept(tmp_path):
    """Pf3D7_07_v3:403266 in a 374-sample sWGA callset: RPBZ 32.2 over 50348 ref and 5886
    alt reads, an effect of 0.128 -- a real pfcrt allele the plain z rule threw out."""
    counts = _strands(21132, 29216, 2555, 3331)
    assert round(_eff(32.2012, 50348, 5886), 3) == 0.128
    v = _vcf(tmp_path, [(1000, _with(RPBZ=32.2012, **counts))])
    a, b, c = (str(tmp_path / f"{x}.bcf") for x in "abc")
    F.hard_qc_filter(v, a, caller="bcftools")
    assert _kept(a) == [1000]
    F.hard_qc_filter(v, b, caller="bcftools", bias_eff=0.1)     # a stricter bar fails it
    assert _kept(b) == []
    F.hard_qc_filter(v, c, caller="bcftools", bias_eff=None)    # and so does z alone
    assert _kept(c) == []


def test_the_effect_gate_never_rescues_a_low_depth_site_the_z_would_fail(tmp_path):
    """Where reads are few, z > 5 already implies a large effect, so the combined rule is
    the z rule: nothing that failed before at single-sample depth passes now."""
    counts = _strands(50, 50, 50, 50)                           # 100 ref, 100 alt
    assert _eff(5.01, 100, 100) > 0.15
    v = _vcf(tmp_path, [(1000, _with(RPBZ=5.01, **counts)), (2000, _with(SCBZ=-5.01, **counts)),
                        (3000, _with(RPBZ=4.99, **counts))])
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(v, out, caller="bcftools")
    assert _kept(out) == [3000]


def test_mapping_quality_bias_is_gated_the_same_way(tmp_path):
    """MQBZ compares ref with alt reads like RPBZ; MQSBZ compares the two strands, so its
    counts are the strand totals."""
    n = 20000
    allele = _strands(n // 2, n // 2, n // 2, n // 2)          # 20000 ref, 20000 alt
    z_small = _z_for(0.05, n, n)                                # ~ 12 at this depth
    assert z_small > 5
    v = _vcf(tmp_path, [(1000, _with(MQBZ=z_small, **allele)),
                        (2000, _with(MQSBZ=z_small, **allele)),
                        (3000, _with(MQBZ=_z_for(0.3, n, n), **allele))])
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(v, out, caller="bcftools")
    assert _kept(out) == [1000, 2000]
    # a lopsided strand table changes MQSBZ's counts but not MQBZ's
    lop = _strands(19000, 1000, 19000, 1000)                    # 38000 fwd, 2000 rev
    z_mqs = _z_for(0.05, 38000, 2000)                           # the same small effect
    v2 = _vcf(tmp_path, [(1000, _with(MQSBZ=z_mqs, **lop)),
                         (2000, _with(MQSBZ=_z_for(0.3, 38000, 2000), **lop))], name="l.vcf")
    out2 = str(tmp_path / "o2.bcf")
    F.hard_qc_filter(v2, out2, caller="bcftools")
    assert _kept(out2) == [1000]


def test_a_non_variant_record_is_not_judged_on_z(tmp_path):
    """No alt reads, no effect to size: the z tests leave a no-ALT record to no_alt_filter,
    the same way the strand-bias test does."""
    hdr = _vcf(tmp_path, [(1000, _with(RPBZ=-14.7, MQBZ=9.0))])
    txt = pathlib.Path(hdr).read_text().replace(
        "chr1\t1000\t.\tA\tT\t222\t.\tDP=40",
        "chr1\t1000\t.\tA\t.\t222\t.\tDP=40").replace("ADF=10,10", "ADF=10").replace(
        "ADR=10,10", "ADR=10")
    v = tmp_path / "noalt.vcf"
    v.write_text(txt)
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(str(v), out, caller="bcftools")
    assert _kept(out) == [1000]
    # bias_eff=None is the old behaviour, and it does remove the record
    out2 = str(tmp_path / "o2.bcf")
    F.hard_qc_filter(str(v), out2, caller="bcftools", bias_eff=None)
    assert _kept(out2) == []


def test_switching_a_z_test_off_switches_its_effect_gate_off_too(tmp_path):
    """read_pos_z=None means no read-position test at all, as it always has -- the effect
    size qualifies a z, it is not a test of its own. With every z test off, ADF/ADR are
    not demanded on its behalf either."""
    v = _vcf(tmp_path, [(1000, _with(RPBZ=40.0, MQBZ=40.0))])
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(v, out, caller="bcftools", read_pos_z=None, max_bias_z=None)
    assert _kept(out) == [1000]
    bare = _vcf(tmp_path, [(1000, {k: val for k, val in CLEAN.items()
                                    if k not in ("ADF", "ADR")})],
                tags=[t for t in BCF_INFO if t not in ("ADF", "ADR")], name="bare.vcf")
    out2 = str(tmp_path / "o2.bcf")
    F.hard_qc_filter(bare, out2, caller="bcftools", sor=None, read_pos_z=None,
                     max_bias_z=None)
    assert _kept(out2) == [1000]


def test_the_effect_gate_names_the_strand_tags_it_reads(tmp_path):
    bare = _vcf(tmp_path, [(1000, {k: val for k, val in CLEAN.items()
                                    if k not in ("ADF", "ADR")})],
                tags=[t for t in BCF_INFO if t not in ("ADF", "ADR")])
    with pytest.raises(SystemExit) as e:
        F.hard_qc_filter(bare, str(tmp_path / "o.bcf"), caller="bcftools", sor=None)
    assert "--bias-eff" in str(e.value)


# --- multiallelic records: the counts behind SOR and the *BZ effect sizes ------------

def _multi_vcf(tmp_path, records, name="multi.vcf"):
    """Like `_vcf`, but each record carries its own ALT column so it can be multiallelic."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>"]
    for t in BCF_INFO:
        n, ty = BCF_INFO[t]
        hdr.append(f'##INFO=<ID={t},Number={n},Type={ty},Description="{t}">')
    hdr.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
    hdr.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1")
    for pos, alt, info in records:
        kv = ";".join(f"{k}={v}" for k, v in info.items())
        hdr.append(f"chr1\t{pos}\t.\tA\t{alt}\t222\t.\t{kv}\tGT\t0/1")
    p = tmp_path / name
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def test_the_bz_effect_size_counts_every_alt_not_just_the_first(tmp_path):
    """Two records with the same ref depth and the same *total* non-reference depth.

    The only difference is that one splits its non-reference reads across two alternates.
    Counting `n2` as ALT1 alone made that record's effect size look larger than it is --
    `eff = z*sqrt((n1+n2+1)/(12*n1*n2))` grows as n2 shrinks -- so the `-e` expression fired
    on it and not on its biallelic twin. The bias runs in the direction that **deletes**
    multiallelic sites, which is the opposite of what this work is for.
    """
    src = _multi_vcf(tmp_path, [
        # ref 50/50, ALT1 10/10, ALT2 90/90 -> non-ref depth 200
        (100, "C,G", {"ADF": "50,10,90", "ADR": "50,10,90", "RPBZ": "4.0",
                      "DP": "400", "MQ": "60", "FS": "0.5", "SCBZ": "0", "MQBZ": "0",
                      "MQSBZ": "0", "BQBZ": "0", "MQ0F": "0"}),
        # ref 50/50, ALT 100/100 -> the same non-ref depth 200, in one allele
        (200, "C", {"ADF": "50,100", "ADR": "50,100", "RPBZ": "4.0",
                    "DP": "400", "MQ": "60", "FS": "0.5", "SCBZ": "0", "MQBZ": "0",
                    "MQSBZ": "0", "BQBZ": "0", "MQ0F": "0"}),
    ])
    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(src, out, caller="bcftools", read_pos_z=3.0, bias_eff=0.15)
    assert _kept(out) == [100, 200], \
        "the multiallelic record must not be judged on one allele's reads"


def test_the_strand_counts_behind_mqsbz_cover_every_allele(tmp_path):
    """MQSBZ compares strands, so its two groups are all-forward and all-reverse reads."""
    src = _multi_vcf(tmp_path, [
        (100, "C,G", {"ADF": "50,10,90", "ADR": "50,10,90", "MQSBZ": "4.0",
                      "DP": "400", "MQ": "60", "FS": "0.5", "RPBZ": "0", "SCBZ": "0",
                      "MQBZ": "0", "BQBZ": "0", "MQ0F": "0"}),
        (200, "C", {"ADF": "50,100", "ADR": "50,100", "MQSBZ": "4.0",
                    "DP": "400", "MQ": "60", "FS": "0.5", "RPBZ": "0", "SCBZ": "0",
                    "MQBZ": "0", "BQBZ": "0", "MQ0F": "0"}),
    ])
    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(src, out, caller="bcftools", max_bias_z=3.0, bias_eff=0.15)
    assert _kept(out) == [100, 200]


def test_sor_reads_the_pooled_non_reference_strand_counts(tmp_path):
    """A record whose *second* alternate is the strand-biased one must still be caught.

    Before this, SOR was computed from ALT1's ADF/ADR alone, so a clean first alternate
    beside a badly strand-skewed second one passed untested.
    """
    src = _multi_vcf(tmp_path, [
        # ALT1 balanced; ALT2 is 200 forward / 0 reverse -- pooled, that is a hard skew
        (100, "C,G", {"ADF": "100,20,200", "ADR": "100,20,0", "DP": "440", "MQ": "60",
                      "FS": "0.5", "RPBZ": "0", "SCBZ": "0", "MQBZ": "0", "MQSBZ": "0",
                      "BQBZ": "0", "MQ0F": "0"}),
        # a clean control with the same totals but no skew
        (200, "C,G", {"ADF": "100,20,100", "ADR": "100,20,100", "DP": "440", "MQ": "60",
                      "FS": "0.5", "RPBZ": "0", "SCBZ": "0", "MQBZ": "0", "MQSBZ": "0",
                      "BQBZ": "0", "MQ0F": "0"}),
    ])
    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(src, out, caller="bcftools", sor=3.0)
    assert _kept(out) == [200], "the skew in the second alternate must be seen"


def test_a_biallelic_record_is_judged_exactly_as_before(tmp_path):
    """The pooled counts must reduce to the old ones at a biallelic site, or every
    existing verdict moves."""
    from plasgenomicsutils.lib.vcf_filters import _N_ALT, _N_FWD, _N_REF, _N_REV

    src = _vcf(tmp_path, [(100, {"ADF": "10,90", "ADR": "10,90", "DP": "200", "MQ": "60",
                                 "FS": "0.5", "RPBZ": "0", "SCBZ": "0", "MQBZ": "0",
                                 "MQSBZ": "0", "BQBZ": "0", "MQ0F": "0"})])
    for expr, want in ((_N_REF, 20), (_N_ALT, 180), (_N_FWD, 100), (_N_REV, 100)):
        got = subprocess.run(["bcftools", "query", "-f", "%POS\n", "-i", f"{expr} == {want}",
                              src], stdout=subprocess.PIPE, text=True,
                             stderr=subprocess.PIPE)
        assert got.stdout.split() == ["100"], f"{expr} should be {want}: {got.stderr[:120]}"


# ---- a declared tag nobody carries ------------------------------------------------------

def test_a_tag_declared_but_missing_on_every_record_is_an_error_naming_it(tmp_path):
    """The header check's trap in another form: `MQ < 55` against a file where every MQ is
    `.` tests nothing and keeps everything."""
    rows = [(1000, {k: v for k, v in CLEAN.items() if k != "MQ"}),
            (2000, {k: v for k, v in CLEAN.items() if k != "MQ"})]
    vcf = _vcf(tmp_path, rows)
    with pytest.raises(SystemExit) as e:
        F.hard_qc_filter(vcf, str(tmp_path / "o.bcf"), caller="bcftools")
    assert "INFO/MQ (read by --mq)" in str(e.value) and "none of the 2 record(s)" in str(e.value)
    # switching that threshold off is enough to proceed
    out = str(tmp_path / "ok.bcf")
    F.hard_qc_filter(vcf, out, caller="bcftools", mq=None)
    assert _kept(out) == [1000, 2000]


def test_a_tag_missing_on_some_records_passes_them_untested_and_is_counted(tmp_path, capsys):
    rows = [(1000, CLEAN), (2000, {k: v for k, v in CLEAN.items() if k != "MQ"}),
            (3000, _with(MQ=30))]
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(_vcf(tmp_path, rows), out, caller="bcftools")
    assert _kept(out) == [1000, 2000]                        # 2000 untested on MQ, 3000 failed it
    cap = capsys.readouterr()
    assert "record(s) without a value pass that test untested -- INFO/MQ 1 of 3" in cap.out + cap.err


def test_a_nan_value_counts_as_missing(tmp_path, capsys):
    """GATK writes MQ=NaN where it could not compute one; `MQ < 55` is false against NaN
    exactly as against `.`, so it is counted the same way."""
    rows = [(1000, CLEAN), (2000, _with(MQ="nan"))]
    out = str(tmp_path / "o.bcf")
    F.hard_qc_filter(_vcf(tmp_path, rows), out, caller="bcftools")
    assert _kept(out) == [1000, 2000]
    cap = capsys.readouterr()
    assert "INFO/MQ 1 of 2" in cap.out + cap.err


def test_gatk_mode_qd_default_is_10(tmp_path):
    """QD 20 was mostly a depth filter by proxy on sWGA data; 10 sits in the valley between
    what VQSR passes and fails on clean biallelic sites (docs/qd_threshold.md)."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>"]
    for t in ("QD", "MQ", "SOR", "MQRankSum", "ReadPosRankSum"):
        hdr.append(f'##INFO=<ID={t},Number=1,Type=Float,Description="{t}">')
    hdr.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
    hdr.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1")
    for pos, qd in ((1000, 9.9), (2000, 10.0), (3000, 15.0), (4000, 30.0)):
        hdr.append(f"chr1\t{pos}\t.\tA\tT\t222\t.\tQD={qd};MQ=60;SOR=1;MQRankSum=0;ReadPosRankSum=0\tGT\t0/1")
    v = tmp_path / "g.vcf"
    v.write_text("\n".join(hdr) + "\n")
    out = str(tmp_path / "o.vcf")
    F.hard_qc_filter(str(v), out)
    assert _kept(out) == [2000, 3000, 4000]
