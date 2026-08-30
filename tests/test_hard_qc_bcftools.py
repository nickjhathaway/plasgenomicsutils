"""hard_qc_filter --caller bcftools: the same questions asked of a bcftools callset.

GATK's metrics are not in a bcftools VCF, so the filter reads what `bcftools mpileup`
writes instead. Two of the mappings are not straight renames, and both are covered here:
bcftools FS is the strand-bias p-value rather than a Phred score, and the *BZ tags are
two-sided.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import vcf_filters as F

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"),
                                reason="bcftools not on PATH")

BCF_INFO = {"DP": ("1", "Integer"), "MQ": ("1", "Integer"), "FS": ("1", "Float"),
            "RPBZ": ("1", "Float"), "SCBZ": ("1", "Float"), "MQBZ": ("1", "Float"),
            "MQSBZ": ("1", "Float"), "BQBZ": ("1", "Float"), "MQ0F": ("1", "Float")}


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


CLEAN = {"DP": 40, "MQ": 60, "FS": 0.25, "RPBZ": 0.0, "SCBZ": 0.0, "MQBZ": 0.0,
         "MQSBZ": 0.0, "BQBZ": 0.0, "MQ0F": 0.0}


def _with(**kw):
    d = dict(CLEAN)
    d.update(kw)
    return d


# ---- the two metrics the filter exists for ---------------------------------------

def test_strand_bias_is_a_p_value_so_small_is_bad(tmp_path):
    # bcftools FS is Fisher's p-value, not GATK's Phred score: a biased site has a tiny
    # FS. Filtering it the GATK way round (FS > 60) would keep exactly the wrong records.
    vcf = _vcf(tmp_path, [(1000, CLEAN), (2000, _with(FS=7.2e-12))])
    out = str(tmp_path / "out.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools")
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
    F.hard_qc_filter(vcf, out, caller="bcftools", read_pos_z=None)
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

def test_qd_is_off_by_default_for_bcftools_but_20_for_gatk(tmp_path):
    # QUAL 222 over DP 40 is 5.6: GATK's QD >= 20 would discard a good bcftools callset,
    # so "auto" leaves it off there
    vcf = _vcf(tmp_path, [(1000, CLEAN)])
    out = str(tmp_path / "auto.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools")
    assert _kept(out) == [1000]
    # asking for it explicitly still works, and does drop the record
    out2 = str(tmp_path / "qd.vcf")
    F.hard_qc_filter(vcf, out2, caller="bcftools", qd=20)
    assert _kept(out2) == []
    # and a modest threshold on the right scale keeps it
    out3 = str(tmp_path / "qd2.vcf")
    F.hard_qc_filter(vcf, out3, caller="bcftools", qd=2)
    assert _kept(out3) == [1000]


def test_depth_resolves_to_info_dp(tmp_path):
    """A bare `DP` in a bcftools expression means FORMAT/DP where both exist."""
    assert "QUAL/INFO/DP" in F._bcftools_qc_expr(
        qd=20, mq=None, strand_bias_p=None, max_bias_z=None, read_pos_z=None,
        mq0f=None, bqbz_z=None)


# ---- refusing to run rather than silently keeping everything ----------------------

def test_a_callset_without_the_tags_is_an_error_naming_them(tmp_path):
    # plain `bcftools mpileup | call` emits RPBZ and MQBZ but not FS, and a comparison
    # against an absent tag is false -- so the filter would quietly keep everything
    vcf = _vcf(tmp_path, [(1000, {"DP": 40, "MQ": 60, "RPBZ": 0.0})],
               tags=("DP", "MQ", "RPBZ"))
    with pytest.raises(SystemExit) as e:
        F.hard_qc_filter(vcf, str(tmp_path / "out.vcf"), caller="bcftools")
    msg = str(e.value)
    assert "INFO/FS" in msg and "INFO/SCBZ" in msg
    assert "bcftools mpileup" in msg          # the call that would produce them
    # ...and disabling the thresholds that read the missing tags is enough to proceed
    out = str(tmp_path / "ok.vcf")
    F.hard_qc_filter(vcf, out, caller="bcftools", strand_bias_p=None, read_pos_z=None,
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
                         mq=None, strand_bias_p=None, read_pos_z=None, max_bias_z=None)


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
    rows = [(1000, "T", CLEAN), (2000, "T", _with(FS=1e-12)),
            (3000, ".", CLEAN), (4000, ".", _with(RPBZ=-14.7))]
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
    # 4000 is non-variant AND read-position biased, so QC takes it out along with 2000
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
