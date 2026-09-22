"""GATK's balanced indel pairs, rewritten as the SNPs they encode.

HaplotypeCaller's haplotype-to-reference alignment (match 200, mismatch -150, gap open -260,
extend -11) scores an insertion plus a deletion above three mismatches, so a multi-base
substitution can arrive as two indel records: pfcrt codons 74-76, ATG AAT AAA to ATT GAA ACA,
come out as `403618 A>AT`, `403622 AT>A`, `403625 A>C`. The records share a PID, which is what
makes the rewrite safe and what this step keys on.
"""

import pathlib
import shutil
import subprocess

import pytest

pytest.importorskip("pysam")
from plasgenomicsutils.lib import shifted_indels as SI
from plasgenomicsutils.lib import filter_pipeline as P

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not on PATH")

# codons 74-76 of pfcrt and some flanking bases; position 1 is the first base below
REF = "GATGAATAAAGGTTC"


@pytest.fixture
def ref_fa(tmp_path):
    p = tmp_path / "ref.fa"
    p.write_text(">chr1\n" + REF + "\n")
    import pysam                                   # pysam ships its own faidx
    pysam.faidx(str(p))
    return str(p)


def _vcf(tmp_path, records, samples, name="in.vcf"):
    """records: (pos, ref, alts). samples: {name: [(gt, pgt) per record]}."""
    hdr = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={len(REF)}>",
           '##INFO=<ID=QD,Number=1,Type=Float,Description="QD">',
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
           '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">',
           '##FORMAT=<ID=PGT,Number=1,Type=String,Description="PGT">',
           '##FORMAT=<ID=PID,Number=1,Type=String,Description="PID">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    pid = f"{records[0][0]}_{records[0][1]}_{records[0][2][0]}"
    for i, (pos, ref, alts) in enumerate(records):
        cols = []
        for s in samples:
            gt, pgt = samples[s][i]
            if gt is None:
                cols.append("./.:0,0:0:.:.")
                continue
            n_alt = 0 if gt == "0/0" else 20
            cols.append(f"{gt}:{40 - n_alt},{n_alt}:40:{pgt or '.'}:{pid if pgt else '.'}")
        hdr.append(f"chr1\t{pos}\t.\t{ref}\t{','.join(alts)}\t500\t.\tQD=25\t"
                   f"GT:AD:DP:PGT:PID\t" + "\t".join(cols))
    p = tmp_path / name
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _sites(path):
    out = subprocess.run(["bcftools", "query", "-f", "%POS\t%REF\t%ALT\n", path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return [tuple(l.split("\t")) for l in out.stdout.splitlines()]


def _gts(path, pos):
    out = subprocess.run(["bcftools", "query", "-i", f"POS=={pos}", "-f", "[%GT ]\n", path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return out.stdout.split()


# the pfcrt cluster: insertion at 2, deletion at 6, SNP at 9 -> substitutions at 4,5,7,9
CVIET = [(2, "A", ("AT",)), (6, "AT", ("A",)), (9, "A", ("C",))]


def test_the_pfcrt_block_becomes_four_snps(tmp_path, ref_fa):
    hom = [("1/1", "1|1")] * 3
    het = [("0/1", "0|1")] * 3
    ref = [("0/0", None)] * 3
    v = _vcf(tmp_path, CVIET, {"carrier_hom": hom, "carrier_het": het, "noncarrier": ref})
    out = str(tmp_path / "o.bcf")
    counts = SI.resolve_shifted_indels(v, out, reference=ref_fa)
    assert counts["resolved"] == 1 and counts["records_removed"] == 3
    assert _sites(out) == [("4", "G", "T"), ("5", "A", "G"), ("7", "T", "A"), ("9", "A", "C")]
    # the carriers are the same at every derived SNP, and the non-carrier is reference
    # (derived genotypes are phased: they were read off a haplotype)
    for pos in (4, 5, 7, 9):
        assert _gts(out, pos) == ["1|1", "0|1", "0|0"]


def test_a_cluster_with_a_net_length_change_is_left_alone(tmp_path, ref_fa):
    """A real indel: an insertion with no matching deletion."""
    recs = [(2, "A", ("AT",)), (9, "A", ("C",))]
    v = _vcf(tmp_path, recs, {"s1": [("0/1", "0|1")] * 2})
    out = str(tmp_path / "o.bcf")
    counts = SI.resolve_shifted_indels(v, out, reference=ref_fa)
    assert counts["resolved"] == 0 and counts["skipped_net_indel"] == 1
    assert _sites(out) == [("2", "A", "AT"), ("9", "A", "C")]


def test_a_snp_only_phase_group_is_untouched(tmp_path, ref_fa):
    recs = [(4, "G", ("T",)), (9, "A", ("C",))]
    v = _vcf(tmp_path, recs, {"s1": [("0/1", "0|1")] * 2})
    out = str(tmp_path / "o.bcf")
    counts = SI.resolve_shifted_indels(v, out, reference=ref_fa)
    assert counts["resolved"] == 0
    assert _sites(out) == [("4", "G", "T"), ("9", "A", "C")]


def test_records_that_share_no_pid_are_not_a_cluster(tmp_path, ref_fa):
    """Unphased records at the same positions: no evidence they are one haplotype."""
    v = _vcf(tmp_path, CVIET, {"s1": [("0/1", None)] * 3})
    out = str(tmp_path / "o.bcf")
    assert SI.resolve_shifted_indels(v, out, reference=ref_fa)["resolved"] == 0
    assert len(_sites(out)) == 3


def _unphased_cluster(tmp_path, n_bad):
    """A cluster where `n_bad` samples are heterozygous at two records with PGT missing on
    one, so which alternate sits on which haplotype is unknown for them."""
    good = [("0/1", "0|1")] * 3
    bad = [("0/1", "0|1"), ("0/1", None), ("0/1", "0|1")]
    samples = {"ok": good}
    samples.update({f"unphased{i}": bad for i in range(n_bad)})
    return _vcf(tmp_path, CVIET, samples)


def test_one_unresolvable_sample_is_the_default_trade(tmp_path, ref_fa):
    """The default is 1: the cluster is rewritten and that sample is set missing. The trade
    is very uneven in real data -- one sample buys hundreds of SNPs, while the worst
    clusters cost dozens of samples for a handful -- which is why the bar sits here."""
    out = str(tmp_path / "o.bcf")
    counts = SI.resolve_shifted_indels(_unphased_cluster(tmp_path, 1), out, reference=ref_fa)
    assert counts["resolved"] == 1 and counts["samples_unresolved"] == 1
    assert _gts(out, 4) == ["0|1", "./."]


def test_a_cluster_costing_more_than_the_bar_is_left_alone(tmp_path, ref_fa):
    """Two unresolvable samples is over the default, so the records stay as they were."""
    out = str(tmp_path / "o.bcf")
    counts = SI.resolve_shifted_indels(_unphased_cluster(tmp_path, 2), out, reference=ref_fa)
    assert counts["resolved"] == 0
    assert counts["skipped_unresolved"] == 1 and counts["samples_kept_by_skipping"] == 2
    assert _sites(out) == [("2", "A", "AT"), ("6", "AT", "A"), ("9", "A", "C")]


def test_zero_never_writes_a_missing_genotype(tmp_path, ref_fa):
    """The strict setting: one unresolvable sample is enough to leave the cluster alone."""
    out = str(tmp_path / "o.bcf")
    counts = SI.resolve_shifted_indels(_unphased_cluster(tmp_path, 1), out, reference=ref_fa,
                                       max_unresolved=0)
    assert counts["resolved"] == 0 and counts["samples_kept_by_skipping"] == 1
    assert _sites(out) == [("2", "A", "AT"), ("6", "AT", "A"), ("9", "A", "C")]


def test_a_sample_with_no_call_at_all_costs_nothing(tmp_path, ref_fa):
    """An unresolvable sample that had no genotype to begin with loses nothing, so the
    cluster is rewritten even at the default."""
    good = [("1/1", "1|1")] * 3
    absent = [(None, None)] * 3
    v = _vcf(tmp_path, CVIET, {"ok": good, "nocall": absent})
    out = str(tmp_path / "o.bcf")
    counts = SI.resolve_shifted_indels(v, out, reference=ref_fa)
    assert counts["resolved"] == 1 and counts["samples_unresolved"] == 0
    assert _gts(out, 4) == ["1|1", "./."]


def test_the_derived_records_carry_depth_info_and_a_provenance_tag(tmp_path, ref_fa):
    v = _vcf(tmp_path, CVIET, {"s1": [("1/1", "1|1")] * 3})
    out = str(tmp_path / "o.bcf")
    SI.resolve_shifted_indels(v, out, reference=ref_fa)
    q = subprocess.run(["bcftools", "query", "-f", "%POS\t%QD\t%SHIFTED\t[%AD:%DP]\n", out],
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout
    for line in q.splitlines():
        pos, qd, shifted, fmt = line.split("\t")
        assert qd == "25" and shifted == "2,6,9" and fmt == "20,20:40"
    # the likelihoods of records that no longer exist are gone
    hdr = subprocess.run(["bcftools", "view", "-h", out], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True).stdout
    assert "ID=PID" not in hdr and "ID=PGT" not in hdr


def test_a_cluster_wider_than_max_span_is_left_alone(tmp_path, ref_fa):
    v = _vcf(tmp_path, CVIET, {"s1": [("0/1", "0|1")] * 3})
    out = str(tmp_path / "o.bcf")
    assert SI.resolve_shifted_indels(v, out, reference=ref_fa, max_span=3)["skipped_span"] == 1
    assert len(_sites(out)) == 3


def test_a_callset_without_phasing_is_copied_through_with_a_note(tmp_path, ref_fa, capsys):
    hdr = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={len(REF)}>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
           "chr1\t4\t.\tG\tT\t500\t.\t.\tGT\t0/1"]
    v = tmp_path / "nophase.vcf"
    v.write_text("\n".join(hdr) + "\n")
    out = str(tmp_path / "o.bcf")
    assert SI.resolve_shifted_indels(str(v), out, reference=ref_fa)["resolved"] == 0
    assert len(_sites(out)) == 1
    cap = capsys.readouterr()
    assert "no FORMAT/PID in the header" in cap.out + cap.err


def test_a_missing_reference_is_an_error_naming_the_flag(tmp_path):
    v = _vcf(tmp_path, CVIET, {"s1": [("0/1", "0|1")] * 3})
    with pytest.raises(SystemExit, match="reference FASTA is needed"):
        SI.resolve_shifted_indels(v, str(tmp_path / "o.bcf"))


def test_the_step_runs_first_in_the_default_chain():
    """Position matters: a cluster only balances while all of its records are present, and
    every later step removes some. On a 249-sample cohort, first in the chain it resolved
    491 clusters into 3,032 SNPs; after the QC and repeat steps, 47 into 251."""
    names = [s["name"] for s in P.DEFAULT_CONFIG["steps"]]
    assert names.index("resolve_shifted_indels") == names.index("no_alt_filter") + 1
    for later in ("hard_qc_filter", "tandem_repeat_mask", "singleton_filter_add_ads",
                  "biallelic_snp_filter"):
        assert names.index("resolve_shifted_indels") < names.index(later)
    assert "resolve_shifted_indels" in P.STEPS
    P.validate_config(P.DEFAULT_CONFIG)


# ---- failing in the first second, not at step nine -------------------------------------

def _gatk_shaped(tmp_path, name="p.vcf", reference=None):
    """A callset with FORMAT/PID, so the step has work to do."""
    v = _vcf(tmp_path, CVIET, {"s1": [("0/1", "0|1")] * 3}, name=name)
    if reference is not None:
        txt = pathlib.Path(v).read_text().replace(
            "##fileformat=VCFv4.2\n", f"##fileformat=VCFv4.2\n##reference={reference}\n", 1)
        pathlib.Path(v).write_text(txt)
    return v


def test_the_pipeline_refuses_up_front_when_the_step_has_no_reference(tmp_path):
    """The failure used to come at the step, after every earlier step had written output."""
    v = _gatk_shaped(tmp_path)
    cfg = {"steps": [{"name": "no_alt_filter"},
                     {"name": "resolve_shifted_indels", "params": {"reference": None}}]}
    with pytest.raises(SystemExit) as e:
        P.run_pipeline(v, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    msg = str(e.value)
    assert "resolve_shifted_indels" in msg and "has no reference" in msg
    assert '"enabled": false' in msg                      # and how to proceed without it
    assert not (tmp_path / "run" / "01_no_alt_filter.bcf").exists()   # nothing was written


def test_a_reference_that_does_not_exist_is_named(tmp_path):
    v = _gatk_shaped(tmp_path)
    cfg = {"steps": [{"name": "resolve_shifted_indels",
                      "params": {"reference": "/no/such/ref.fa"}}]}
    with pytest.raises(SystemExit, match="does not exist: /no/such/ref.fa"):
        P.run_pipeline(v, str(tmp_path / "run"), cfg, emit_snp_bed=False)


def test_a_header_reference_that_is_not_on_this_machine_is_named(tmp_path):
    """A callset called elsewhere carries its own ##reference path, which is the usual way
    this bites: the config looks fine and the path is somebody else's."""
    v = _gatk_shaped(tmp_path, reference="/elsewhere/Pf3D7.fasta")
    cfg = {"steps": [{"name": "resolve_shifted_indels"}]}
    with pytest.raises(SystemExit, match="not on this machine"):
        P.run_pipeline(v, str(tmp_path / "run"), cfg, emit_snp_bed=False)


def test_the_headers_reference_is_accepted_when_it_is_there(tmp_path, ref_fa):
    v = _gatk_shaped(tmp_path, reference=ref_fa)
    cfg = {"steps": [{"name": "resolve_shifted_indels"}]}
    tally = P.run_pipeline(v, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    assert tally[-1]["variants"] == 4                     # the block was rebuilt


def test_a_callset_with_no_phasing_needs_no_reference(tmp_path):
    """No FORMAT/PID means no clusters to rebuild, so the step is a no-op -- refusing it
    for want of a reference would turn a good bcftools run into an error."""
    hdr = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={len(REF)}>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
           "chr1\t4\t.\tG\tT\t500\t.\t.\tGT\t0/1"]
    v = tmp_path / "nophase.vcf"
    v.write_text("\n".join(hdr) + "\n")
    cfg = {"steps": [{"name": "resolve_shifted_indels"}]}
    tally = P.run_pipeline(str(v), str(tmp_path / "run"), cfg, emit_snp_bed=False)
    assert tally[-1]["variants"] == 1


def test_a_disabled_step_is_not_preflighted(tmp_path):
    v = _gatk_shaped(tmp_path)
    cfg = {"steps": [{"name": "resolve_shifted_indels", "enabled": False},
                     {"name": "no_alt_filter"}]}
    tally = P.run_pipeline(v, str(tmp_path / "run"), cfg, emit_snp_bed=False)
    assert any(t.get("skipped") for t in tally)


def test_the_provenance_tag_is_declared_as_a_list(tmp_path, ref_fa):
    """One position per source record, and a cluster has at least two, so `Number=1` makes
    every derived record malformed. A strict reader refuses it: SeqArray stops the GDS
    build with "INFO ID 'SHIFTED' should have 1 value(s) but receives 2"."""
    out = str(tmp_path / "o.bcf")
    SI.resolve_shifted_indels(_vcf(tmp_path, CVIET, {"s1": [("1/1", "1|1")] * 3}), out,
                              reference=ref_fa)
    hdr = subprocess.run(["bcftools", "view", "-h", out], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True).stdout
    line = next(l for l in hdr.splitlines() if "ID=SHIFTED" in l)
    assert "Number=.," in line and "Number=1," not in line
    # and the value really is several, so the declaration is the honest one
    got = subprocess.run(["bcftools", "query", "-f", "%SHIFTED\n", out],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout
    assert all(v == "2,6,9" for v in got.split())
