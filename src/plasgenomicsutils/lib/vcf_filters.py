"""Parameterized VCF/BCF filtering steps, each backed by bcftools/bedtools.

Every step takes an input and output path plus tunable thresholds, so the same
building blocks can be run individually or chained by ``filter_pipeline``.
"""

from __future__ import annotations

import contextlib
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

from ..utils.small_utils import Utils
from .reporting import detail, say
from .bcftools import (count_variants, count_variants_matching, format_tags, index_vcf,
                       info_tags, out_flag, q, require, sh)
from .strip_format import GENOTYPE_LINKED_FORMAT, strip_stale_format


#: What `bcftools mpileup` has to be asked for so a callset carries everything the
#: ``caller="bcftools"`` filter reads. The rest of the tags it uses (RPBZ, MQBZ, MQSBZ,
#: BQBZ, SCBZ, MQ0F, MQ, DP) are emitted by default and cannot be requested explicitly.
BCFTOOLS_MPILEUP_ANNOTATIONS = (
    "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,FORMAT/SCR,"
    "INFO/AD,INFO/ADF,INFO/ADR,INFO/FS,INFO/SCR"
)

#: What `bcftools mpileup` writes whether asked to or not -- these cannot be requested
#: explicitly (it rejects them in -a), so a check that the filter's tags are obtainable has
#: to know about them. MQ and DP are always present.
BCFTOOLS_MPILEUP_DEFAULT_TAGS = frozenset({
    "BQBZ", "IDV", "IMF", "MQ0F", "MQBZ", "MQSBZ", "RPBZ", "SCBZ", "SGB", "VDB",
    "MQ", "DP",
})

#: The call that produces such a callset, quoted back at anyone whose input lacks the tags.
BCFTOOLS_CALL_RECIPE = (
    f"bcftools mpileup -f REF.fa -a {BCFTOOLS_MPILEUP_ANNOTATIONS} IN.bam -Ou \\\n"
    "    | bcftools call -m -Ob -o OUT.bcf        # or: plasgenomicsutils call_variants"
)


def _sor_from_ad_expr(threshold: float) -> str:
    """``SOR > threshold`` written out of INFO/ADF and INFO/ADR, since bcftools writes no SOR.

    GATK's StrandOddsRatio on the 2x2 table of ref/alt by forward/reverse, with the +1
    pseudocount, is::

        R    = (refFwd/refRev) * (altRev/altFwd)
        SOR  = ln(R + 1/R) + ln(min/max of the ref pair) - ln(min/max of the alt pair)

    which is an **effect size**: scale every count by the same factor and it does not move.
    That is the whole reason to compute it rather than read ``FS``. ``FS`` is a p-value, so
    it answers "am I sure there is a skew" -- a question whose answer is always yes once a
    cohort is large enough, because significance grows with the reads pooled across samples
    while the skew itself stays put.

    bcftools has no logarithm, so the comparison is exponentiated and cross-multiplied into
    integer arithmetic, which is exact and cannot overflow at any plausible depth::

        (A^2 + B^2) * min(ref) * max(alt)  >  e^threshold * A * B * max(ref) * min(alt)

    with ``A = refFwd*altRev``, ``B = refRev*altFwd``. ``min``/``max`` of a pair come from
    ``(a+b-|a-b|)`` and ``(a+b+|a-b|)``; their halves cancel between the two sides.

    ``INFO/`` is spelled out on every tag: a bare ``ADF`` is ambiguous where FORMAT/ADF also
    exists, and bcftools refuses the expression rather than guessing. Only the first ALT is
    tested, which is what the biallelic callsets this chain produces contain.
    """
    ref_f, ref_r = "(INFO/ADF[0]+1)", "(INFO/ADR[0]+1)"
    alt_f, alt_r = "(INFO/ADF[1]+1)", "(INFO/ADR[1]+1)"
    a, b = f"({ref_f}*{alt_r})", f"({ref_r}*{alt_f})"
    lo = lambda x, y: f"({x}+{y}-abs({x}-{y}))"       # noqa: E731 - 2x min(x, y)
    hi = lambda x, y: f"({x}+{y}+abs({x}-{y}))"       # noqa: E731 - 2x max(x, y)
    return (f"(N_ALT>0 && ({a}*{a}+{b}*{b})*{lo(ref_f, ref_r)}*{hi(alt_f, alt_r)} > "
            f"{math.exp(threshold):.10g}*{a}*{b}*{hi(ref_f, ref_r)}*{lo(alt_f, alt_r)})")


#: Allele counts behind the ref-vs-alt z-scores (RPBZ, SCBZ, MQBZ, BQBZ): first ALT only,
#: as for SOR. INFO/AD is these two sums, but ADF/ADR are what SOR already requires.
_N_REF = "(INFO/ADF[0]+INFO/ADR[0])"
_N_ALT = "(INFO/ADF[1]+INFO/ADR[1])"
#: Strand counts behind MQSBZ, which compares mapping quality between strands, not alleles.
_N_FWD = "(INFO/ADF[0]+INFO/ADF[1])"
_N_REV = "(INFO/ADR[0]+INFO/ADR[1])"


def _mwu_bias_expr(tag: str, z: float, eff: float | None, n1: str, n2: str) -> str:
    """``|tag| > z``, and -- unless ``eff`` is None -- the effect behind that z ``> eff``.

    The ``*BZ`` tags are Mann-Whitney z-scores comparing two groups of reads pooled over
    every sample at the site: ref against alt for ``RPBZ``/``SCBZ``/``MQBZ``/``BQBZ``,
    forward against reverse for ``MQSBZ``. A z-score grows with the reads behind it: the
    same modest shift in read position scores z = 1 in one sample and z = 30 in a cohort of
    400, so a fixed z cutoff asks "am I sure there is a shift" -- a question whose answer is
    always yes once a cohort is large enough. The effect size the z is built from does not
    move with depth. With ``n1`` and ``n2`` reads in the two groups::

        eff = z * sqrt((n1 + n2 + 1) / (12 * n1 * n2))

    is how far ``P(a read from one group ranks above a read from the other)`` sits from
    0.5, so ``eff = 0.15`` means a 65:35 split (rank-biserial correlation 0.3). Multiply
    every count by the same factor and it is unchanged -- the property the strand-bias test
    already has through SOR.

    The two are combined rather than swapped because the effect size is noisy where the z
    is not: at a single sample's depth a shift of 0.15 is well within sampling error, and
    a threshold on the effect alone fails four times as many low-depth sites as the z does,
    all of them noise. Requiring both is exactly the z rule where reads are few (z = 5 at
    100 reads per group already implies eff = 0.2) and exactly the effect rule where reads
    are many, with no depth at which the combined rule is stricter than the z alone.

    bcftools has no square root, so the comparison is squared and cross-multiplied::

        z^2 * (n1 + n2 + 1) > 12 * eff^2 * n1 * n2

    A record whose ADF/ADR carry no ALT entry (a no-ALT site) makes the effect clause
    false, so with ``eff`` set these tests do not judge non-variant records;
    ``no_alt_filter`` removes those on its own terms.
    """
    parts = [f'{tag}!="."', f"abs({tag}) > {z}"]
    if eff is not None:
        # N_ALT>0 first: on a no-ALT record ADF[1]/ADR[1] do not exist, and bcftools does
        # not treat arithmetic on a missing index as false -- the same guard SOR carries
        parts.append("N_ALT>0")
        parts.append(f"{tag}*{tag}*({n1}+{n2}+1) > {12 * eff * eff:.10g}*{n1}*{n2}")
    return "(" + " && ".join(parts) + ")"


def _bcftools_qc_expr(*, qd, mq, sor, strand_bias_p, max_bias_z, read_pos_z, mq0f, bqbz_z,
                      bias_eff):
    """The bcftools-native mirror of the GATK hard filter, as a bcftools -e expression.

    Each GATK metric has a counterpart in what `bcftools mpileup` writes, but two of them
    do not carry over naively:

    * GATK ``FS`` is Phred-scaled and rises with bias; **bcftools ``FS`` is the p-value
      itself**, so the test is ``FS <`` a small number rather than ``FS >`` a large one.
      (Phred 60, GATK's usual cutoff, is p = 1e-6 -- the same statement either way.)
      Both are p-values, and both therefore measure the cohort as much as the site, which
      is why ``sor`` and not ``strand_bias_p`` is the strand-bias test on by default --
      see :func:`_sor_from_ad_expr`.
    * GATK's rank-sum filters are one-sided, because the sign says which way the alt
      allele leans. The bcftools ``*BZ`` tags are documented as "closer to 0 is better",
      and a read-position artifact shows up as either sign, so they are tested on
      ``abs()``.
    * Every ``*BZ`` tag is a z-score, and a z-score grows with the reads pooled across
      samples. Each is therefore tested together with the effect size behind it, computed
      from ``INFO/ADF``/``INFO/ADR`` -- allele counts for the ref-vs-alt tags, strand
      counts for ``MQSBZ`` -- unless ``bias_eff`` is None. See :func:`_mwu_bias_expr`.

    ``DP`` is written as ``INFO/DP`` deliberately: a bare ``DP`` in a bcftools expression
    resolves to ``FORMAT/DP`` where both exist, which silently changes what is filtered.
    """
    parts = []
    if qd is not None:
        parts.append(f"QUAL/INFO/DP < {qd}")
    if mq is not None:
        parts.append(f"MQ < {mq}")
    if sor is not None:                                 # strand bias, as an effect size
        parts.append(_sor_from_ad_expr(sor))
    if strand_bias_p is not None:                       # strand bias, Fisher p-value
        parts.append(f'(FS!="." && FS < {strand_bias_p})')
    # each z test is "significant and large" -- see _mwu_bias_expr for why not one or the
    # other; bias_eff=None makes them plain z tests again
    if read_pos_z is not None:                          # the alt sits near the read ends
        parts.append(_mwu_bias_expr("RPBZ", read_pos_z, bias_eff, _N_REF, _N_ALT))
        parts.append(_mwu_bias_expr("SCBZ", read_pos_z, bias_eff, _N_REF, _N_ALT))
    if max_bias_z is not None:                          # mapping quality, and vs strand
        parts.append(_mwu_bias_expr("MQBZ", max_bias_z, bias_eff, _N_REF, _N_ALT))
        parts.append(_mwu_bias_expr("MQSBZ", max_bias_z, bias_eff, _N_FWD, _N_REV))
    if bqbz_z is not None:
        parts.append(_mwu_bias_expr("BQBZ", bqbz_z, bias_eff, _N_REF, _N_ALT))
    if mq0f is not None:
        parts.append(f'(MQ0F!="." && MQ0F > {mq0f})')
    if not parts:
        raise SystemExit("hard_qc_filter: every bcftools threshold was disabled, so this "
                         "step would keep every record; drop it from the chain instead")
    return " || ".join(parts)


#: INFO tags each bcftools-mode threshold reads, so a missing one is named rather than
#: silently never matching -- a bcftools comparison against an absent tag is just false.
_BCFTOOLS_QC_NEEDS = {
    "qd": ("DP",), "strand_bias_p": ("FS",), "read_pos_z": ("RPBZ", "SCBZ"),
    "bias_eff": ("ADF", "ADR"),
    "sor": ("ADF", "ADR"),
    "max_bias_z": ("MQBZ", "MQSBZ"), "bqbz_z": ("BQBZ",), "mq0f": ("MQ0F",),
    "mq": ("MQ",),
}


#: INFO tags each GATK-mode threshold reads.
_GATK_QC_NEEDS = {
    "qd": ("QD",), "mq": ("MQ",), "sor": ("SOR",), "mqranksum": ("MQRankSum",),
    "readposranksum": ("ReadPosRankSum",), "fs": ("FS",),
}


def _check_qc_tags(inp: str, wanted: dict, needs: dict, label: str, advice: str) -> None:
    """Fail before running bcftools when the callset lacks a tag a threshold reads.

    bcftools aborts on an undefined tag, but from inside a shell pipeline: the error names
    the tag and then the run dies with `unknown file type` from the next process in the
    pipe, which reads as a corrupt input rather than a missing annotation. Checking the
    header first turns that into one sentence naming every missing tag at once, and which
    threshold to switch off if the annotation is genuinely unavailable.
    """
    have = info_tags(inp)
    missing = defaultdict(set)
    for k, v in wanted.items():
        if v is None:
            continue
        for t in needs[k]:
            if t not in have:
                missing[t].add("--" + k.replace("_", "-"))
    if not missing:
        return
    # naming the threshold beside the tag, since which flag to reach for is the next
    # question and the mapping from tag to flag is not obvious from either end
    named = ", ".join(f"INFO/{t} (read by {', '.join(sorted(f))})"
                      for t, f in sorted(missing.items()))
    raise SystemExit(
        f"{label}: this callset has no " + named + ".\n"
        "  A comparison against a tag that is not there is simply false, so the filter "
        "would keep\n  everything and say nothing. " + advice
    )


def _check_bcftools_qc_tags(inp: str, wanted: dict) -> None:
    _check_qc_tags(
        inp, wanted, _BCFTOOLS_QC_NEEDS, "hard_qc_filter --caller bcftools",
        "Call with the annotations it needs:\n\n"
        f"    {BCFTOOLS_CALL_RECIPE}\n\n"
        "  ...or disable the thresholds that read the missing tags.")


def _check_gatk_qc_tags(inp: str, wanted: dict) -> None:
    _check_qc_tags(
        inp, wanted, _GATK_QC_NEEDS, "hard_qc_filter --caller gatk",
        "These are GATK's annotations;\n"
        "  a callset from `bcftools call` does not carry them and wants "
        "`--caller bcftools` instead.\n"
        "  A callset that has been through an earlier filter may simply have had its INFO "
        "stripped.\n"
        "  Otherwise set the thresholds that read the missing tags to null to switch them "
        "off.")


def _trim_alt_alleles(inp: str, out: str) -> int:
    """Drop ALT alleles no genotype carries; return how many records lost one.

    A joint callset subset to fewer samples keeps every ALT the full cohort had, now
    carried by nobody: `bcftools view -S` removes samples, not alleles. Those alleles are
    still counted, still classified, and still filtered on, so the numbers describe a
    cohort that is not the one in the file.

    A genotype-linked Number=G field (e.g. PL) whose length disagrees with the genotypes
    makes `--trim-alt-alleles` abort, so the same guard `biallelic_snp_filter` uses applies
    here: null just the inconsistent records first, keeping the valid likelihoods.
    """
    src, tmp = inp, None
    try:
        if any(t in format_tags(inp) for t in GENOTYPE_LINKED_FORMAT):
            tmp = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
            strip_stale_format(inp, tmp, fields=GENOTYPE_LINKED_FORMAT, mode="mismatch")
            src = tmp
        sh(f"bcftools view --trim-alt-alleles {q(src)} -Ob -o {q(out)}", tools=("bcftools",))
    finally:
        if tmp:
            for suffix in ("", ".csi"):
                try:
                    os.remove(tmp + suffix)
                except OSError:
                    pass
    return _count_alt_changes(inp, out)


def _count_alt_changes(before: str, after: str) -> int:
    """How many records' ALT columns differ between two callsets.

    Trimming never drops a record, so the two stream in lockstep and this needs no index
    and no memory beyond a line at a time.
    """
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%ALT\n"]
    a = subprocess.Popen(cmd + [str(before)], stdout=subprocess.PIPE, text=True)
    b = subprocess.Popen(cmd + [str(after)], stdout=subprocess.PIPE, text=True)
    try:
        return sum(1 for x, y in zip(a.stdout, b.stdout) if x != y)
    finally:
        for proc in (a, b):
            proc.stdout.close()
            proc.wait()


def no_alt_filter(inp: str, out: str, *, keep: bool = False, trim: bool = False,
                  keep_bed: str | None = None) -> int:
    """Drop records with no ALT allele, i.e. positions that turned out non-variant.

    Calling a list of positions (``bcftools call`` without ``-v``) reports every one of
    them, so a callset over a region list carries a record wherever a sample is simply
    reference. Removing them in their own step rather than letting a QC rule do it keeps
    the two reasons apart in ``variant_counts.tsv``: how many positions had nothing to
    call, and how many real variants failed quality.

    It matters because **the bias statistics are computed whether or not an ALT was
    called**. A non-variant record still carries FS, RPBZ, MQBZ and the rest, describing
    the non-reference reads that were there but not called -- so a hard QC filter does
    remove non-variant records, and on a real region list it removed a quarter of them.
    Those are positions where the non-reference evidence is strand- or end-biased enough
    that no ALT was called, which is arguably where a reference call is least safe; either
    way, deciding it here rather than as a side effect is the point.

    ``keep`` passes them through, for a fill-in workflow where "this sample is reference
    here" is the answer being sought. The count is reported either way.

    ``trim`` first drops ALT alleles no genotype carries, which is what makes every count
    after this step describe the cohort in the file rather than the one it was subset from.
    The two are reported separately -- how many records lost an allele, and how many were
    left with none -- because they are different facts, and trimming is what creates some
    of the second kind.
    """
    fmt = out_flag(out)
    trimmed = None
    try:
        if trim:
            # trim BEFORE the count: a record whose every ALT was uncarried becomes ALT="."
            # here, and it is this step's job to say so rather than the next step's to
            # wonder where it came from
            trimmed = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
            n_trimmed = _trim_alt_alleles(inp, trimmed)
            inp = trimmed
            sys.stderr.write(f"NOTE: {n_trimmed} record(s) had ALT allele(s) no genotype "
                             f"carries trimmed away\n")
        return _no_alt_filter(inp, out, fmt, keep=keep, keep_bed=keep_bed)
    finally:
        if trimmed and os.path.exists(trimmed):
            os.unlink(trimmed)


def _no_alt_filter(inp: str, out: str, fmt: str, *, keep: bool, keep_bed: str | None) -> int:
    n_no_alt = count_variants_matching(inp, 'ALT="."')
    verb = "kept" if keep else "dropped"
    sys.stderr.write(f"NOTE: {n_no_alt} record(s) have no ALT allele ({verb})\n")
    if keep:
        sh(f"bcftools view {q(inp)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return 0
    expr = 'ALT="."'
    if not keep_bed:
        sh(f"bcftools view -e {q(expr)} {q(inp)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return 0
    prep = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
    try:
        shutil.copyfile(inp, prep)
        sh(f"bcftools view -e {q(expr)} {q(prep)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return _rescue_whitelisted(prep, out, keep_bed, "no_alt_filter")
    finally:
        if os.path.exists(prep):
            os.unlink(prep)


def hard_qc_filter(inp: str, out: str, *, caller: str = "gatk", qd: float | None | str = "auto",
                   mq: float | None = 55, sor: float | None = 3, mqranksum: float = -5.0,
                   readposranksum: float = -5.0, fs: float | None = None,
                   strand_bias_p: float | None | str = "auto", read_pos_z: float | None = 5.0,
                   max_bias_z: float | None = 5.0, bqbz_z: float | None = None,
                   bias_eff: float | None = 0.15,
                   mq0f: float | None = None, keep_bed: str | None = None) -> int:
    """Hard filter on INFO metrics; keep records still flagged PASS.

    ``caller="gatk"`` (the default) filters on QD / MQ / SOR / MQRankSum / ReadPosRankSum,
    which is what GATK writes. ``caller="bcftools"`` filters on the counterparts
    `bcftools mpileup` writes -- see :func:`_bcftools_qc_expr` for what maps to what, and
    where the mapping is not a straight rename.

    ``qd="auto"`` resolves per caller: 20 for GATK, and **off** for bcftools, whose QUAL is
    not on the same scale -- a 40x site called at QUAL 222 has QUAL/DP of 5.6, so carrying
    GATK's 20 across would discard a perfectly good callset. Pass a number to set it
    anyway, or ``None`` to switch it off.

    ``sor`` is the strand-bias test in **both** modes and means the same thing in both: GATK
    writes SOR, and for a bcftools callset it is computed from ``INFO/ADF``/``INFO/ADR`` (see
    :func:`_sor_from_ad_expr`), so the threshold carries across unchanged. ``strand_bias_p``
    tests bcftools' ``FS`` instead and is ``"auto"`` -- meaning off -- because a fixed
    p-value cutoff gets stricter as a cohort grows while the skew it is meant to detect does
    not. Set it to a number to add it back.

    ``bias_eff`` qualifies every z test in bcftools mode (``read_pos_z`` on RPBZ/SCBZ,
    ``max_bias_z`` on MQBZ/MQSBZ, ``bqbz_z`` on BQBZ): a site fails a z test only when the
    z exceeds its cutoff **and** the effect size behind that z exceeds ``bias_eff`` (see
    :func:`_mwu_bias_expr` for the arithmetic and why both). The z alone tightens as a
    cohort grows -- a shift that scores z = 1 in one sample scores z = 30 pooled over 400 --
    and the effect alone is noise at a single sample's depth. ``bias_eff=None`` restores
    the plain z tests. Setting a z cutoff to None switches that test off, as before.

    The GATK-only thresholds (``sor``, ``mqranksum``, ``readposranksum``, ``fs``) are
    ignored in bcftools mode, and the bcftools-only ones in GATK mode; each set names tags
    the other caller does not write.

    ``keep_bed`` whitelists regions from this rule; a rescued record keeps its ``FAIL``
    FILTER, so a variant kept despite failing QC still says that it failed.
    """
    if caller not in ("gatk", "bcftools"):
        raise SystemExit(f"hard_qc_filter: caller must be gatk or bcftools, not {caller!r}")
    # "auto" means "whatever suits this caller": QD 20 is GATK's, and QUAL/DP on bcftools
    # output is not the same quantity, so it is left off there rather than reused.
    if qd == "auto":
        qd = 20.0 if caller == "gatk" else None
    # FS is off by default for the reason given in _sor_from_ad_expr: it is a p-value, so
    # a fixed cutoff tightens as a cohort grows, and `sor` asks the same question of the
    # skew itself. Pass a number to switch it back on.
    if strand_bias_p == "auto":
        strand_bias_p = None

    if caller == "bcftools":
        # the effect size only reads ADF/ADR on behalf of a z test, so with every z test
        # off it must not ask for them
        if read_pos_z is None and max_bias_z is None and bqbz_z is None:
            bias_eff = None
        thresholds = dict(qd=qd, mq=mq, sor=sor, strand_bias_p=strand_bias_p,
                          read_pos_z=read_pos_z, max_bias_z=max_bias_z, bqbz_z=bqbz_z,
                          mq0f=mq0f, bias_eff=bias_eff)
        _check_bcftools_qc_tags(inp, thresholds)
        expr = _bcftools_qc_expr(**thresholds)
    else:
        thresholds = dict(qd=qd, mq=mq, sor=sor, mqranksum=mqranksum,
                          readposranksum=readposranksum, fs=fs)
        _check_gatk_qc_tags(inp, thresholds)
        # each threshold is skipped when it is None, so switching one off is a way to run
        # against a callset that lacks its tag rather than an error with no way forward
        forms = {
            "qd": "QD < {}", "mq": "MQ < {}", "sor": "SOR > {}",
            "mqranksum": '(MQRankSum!="." && MQRankSum < {})',
            "readposranksum": '(ReadPosRankSum!="." && ReadPosRankSum < {})',
            "fs": "FS > {}",
        }
        parts = [forms[k].format(v) for k, v in thresholds.items() if v is not None]
        if not parts:
            raise SystemExit(
                "hard_qc_filter --caller gatk: every threshold is off, so this step would "
                "copy its input; drop it from the chain instead")
        expr = " || ".join(parts)
    fmt = out_flag(out)
    if not keep_bed:
        sh(f"bcftools filter -m + -s FAIL -e {q(expr)} {q(inp)} -Ou "
           f"| bcftools view -f PASS -O{fmt} -o {q(out)}", tools=("bcftools",))
        return 0
    # flagged first, selected second, so the whitelist has a header-compatible source
    prep = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
    try:
        sh(f"bcftools filter -m + -s FAIL -e {q(expr)} {q(inp)} -Ob -o {q(prep)}",
           tools=("bcftools",))
        sh(f"bcftools view -f PASS {q(prep)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return _rescue_whitelisted(prep, out, keep_bed, "hard_qc_filter")
    finally:
        if os.path.exists(prep):
            os.unlink(prep)


def singleton_add_ads(inp: str, out: str, *, min_samples: int = 1,
                     keep_bed: str | None = None) -> int:
    """Drop variants seen as ALT in <= min_samples samples; add FORMAT/ADS.

    ADS (summed allelic depth, the reads actually used to genotype) is a more
    honest per-site depth than DP. It is the sum over *all* AD entries via
    ``smpl_sum``, so it stays correct on multiallelic sites — the older
    ``AD[*:0]+AD[*:1]`` form counted only ref + first ALT and silently
    undercounted anything with a second ALT allele.
    """
    fmt = out_flag(out)
    # Single '&' is required inside COUNT(): it combines the conditions per
    # sample, so COUNT() tallies samples meeting both. '&&' would collapse them
    # across the whole record and match every site.
    keep = f'COUNT(GT!="RR" & GT!="mis") > {min_samples}'
    ads = "FORMAT/ADS=int(smpl_sum(FORMAT/AD))"
    if not keep_bed:
        sh(f"bcftools view -i {q(keep)} {q(inp)} -Ou "
           f"| bcftools +fill-tags -O{fmt} -o {q(out)} -- -t {q(ads)}", tools=("bcftools",))
        return 0
    # ADS is per-sample and does not depend on which records survive, so tagging everything
    # first and selecting second gives the same output and leaves a source to rescue from
    prep = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
    try:
        sh(f"bcftools +fill-tags {q(inp)} -Ob -o {q(prep)} -- -t {q(ads)}",
           tools=("bcftools",))
        sh(f"bcftools view -i {q(keep)} {q(prep)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return _rescue_whitelisted(prep, out, keep_bed, "singleton_filter_add_ads")
    finally:
        if os.path.exists(prep):
            os.unlink(prep)


#: Every variant type that is not a SNP. Selecting SNPs by exclusion rather than with
#: ``-v snps`` is deliberate: ``-v snps`` keeps a record if **any** of its alleles is a SNP,
#: so a mixed ``A>T,ATT`` site passes it. Excluding the rest leaves only records whose every
#: allele is a substitution. ``ref`` is in the list because a record left with no ALT is not
#: a variant, and would otherwise survive when the biallelic test is off.
NON_SNP_TYPES = "indels,mnps,other,ref,bnd"


MNP_HANDLING = ("split", "remove", "keep")


def _snp_select_args(snps_only: bool, biallelic: bool, mnp_handling: str = "remove") -> str:
    parts = []
    if snps_only:
        # `remove` is the only mode that excludes MNPs here. `keep` wants them, and `split`
        # has to let them through so there is something left for the atomiser to split --
        # after which there are none, because each has become its component SNPs.
        types = (NON_SNP_TYPES if mnp_handling == "remove"
                 else ",".join(t for t in NON_SNP_TYPES.split(",") if t != "mnps"))
        parts.append(f"-V {types}")
    # A record with no ALT has to go when the type test is on, and `-V ref` cannot be
    # trusted to do it: whether ALT="." counts as type `ref` changed between bcftools
    # releases. On 1.19 `-v ref` selects nothing at all, so `-V ref` excludes nothing and a
    # non-variant record survives a SNPs-only filter; on 1.24 it is excluded as intended.
    # `--min-alleles 2` asks the same question by counting alleles, which every version
    # answers the same way. `ref` stays in the exclusion list as the statement of intent.
    if snps_only or biallelic:
        parts.append("-m2")
        # A record whose only ALT is `*` states that an upstream deletion covers this
        # position and nothing else -- two alleles by count, no variant to call. It passes
        # every test above, so it needs naming. The N_ALT test is what keeps `A > *,T`, a
        # real SNP that merely sits under a deletion in some samples.
        parts.append("""-e 'N_ALT=1 && ALT="*"'""")
    if biallelic:
        parts.append("-M2")
    return " ".join(parts)


def biallelic_snp_filter(inp: str, out: str, *, trim: bool = True,
                         snps_only: bool = True, biallelic: bool = True,
                         mnp_handling: str = "split") -> None:
    """Keep SNPs and/or biallelic sites, optionally trimming unused ALT alleles first.

    The two questions are separate knobs because they are separate questions. ``snps_only``
    drops anything that is not a substitution; ``biallelic`` drops anything with more than one
    ALT. Both default on, which is the historical behaviour and what most population-genetic
    tooling wants. Turning ``biallelic`` off keeps multiallelic SNPs, for downstream tools
    that can read them -- the site set is otherwise identical.

    ``trim`` (default) drops ALT alleles absent from every genotype before either test
    (``bcftools view --trim-alt-alleles``). Run this AFTER re-genotyping
    (:func:`~plasgenomicsutils.lib.regenotype.filter_ad_regenotype`): a site that looked
    multiallelic only because one sample carried a low-level artifact allele -- now
    re-genotyped away -- collapses to a genuine biallelic SNP and is kept, instead of being
    discarded by a naive ``-m2 -M2``. Trimming is also what turns a mixed ``A>T,ATT`` record
    into a plain SNP when nothing carried the indel, so the type test belongs after it too.

    ``mnp_handling`` decides what happens to a multi-base substitution:

    * ``"split"`` (default) breaks it into its component SNPs with ``bcftools norm -a``,
      which needs no reference. This also rewrites a substitution that was merely *written*
      with padding -- ``REF=TTATA ALT=CTATA`` differs at one base and is a SNP -- into its
      minimal form, which is what makes a callset safe for tools that assume one base per
      SNP record.
    * ``"remove"`` drops the record, the behaviour of the type test alone.
    * ``"keep"`` leaves it as it is, for a downstream tool that reads MNPs.

    ``"split"`` composes with ``biallelic=False``: multiallelic sites are split apart,
    atomised, and rejoined, which leaves them exactly as they were while still breaking up
    any MNP among them.

    With both tests off the step only trims, which it says rather than looking like a no-op.
    """
    if mnp_handling not in MNP_HANDLING:
        raise SystemExit(f"biallelic_snp_filter: mnp_handling must be one of "
                         f"{', '.join(MNP_HANDLING)}, not {mnp_handling!r}")
    fmt = out_flag(out)
    select = _snp_select_args(snps_only, biallelic, mnp_handling)
    split = mnp_handling == "split"
    if not select:
        say("     biallelic_snp_filter: no type or allele-count test is on, so this step "
              + ("only splits multi-base substitutions and trims unused ALT alleles"
                 if split else "only trims unused ALT alleles"))
        if not trim and not split:
            raise SystemExit(
                "biallelic_snp_filter: trim, snps_only and biallelic are all off and "
                "mnp_handling does not split, so this step would copy its input; drop it "
                "from the chain instead")
    # Atomise *after* the selection, never before: the selection decides which records
    # exist, and the atomiser is happy to take a multiallelic record apart as well as an
    # MNP -- one output record per ALT, `*` standing in for the others.
    #
    # Splitting the site first is what stops that. `-m -any` leaves the atomiser nothing
    # but biallelic records, where it can only do what was asked; `-m +any` then puts the
    # sites back together. On a real callset the round trip is exact -- record for record,
    # including INFO and every FORMAT field -- because nothing was ever decomposed for the
    # wrong reason, only re-expressed.
    atomize = (f" | bcftools norm -m -any - -Ou | bcftools norm -a - -Ou"
               f" | bcftools norm -m +any - -O{fmt} -o {q(out)}") if split else ""
    tail = f" -O{fmt} -o {q(out)}" if not split else " -Ou"
    if not trim:
        sh(f"bcftools view {select} {q(inp)}{tail}{atomize}", tools=("bcftools",))
        return
    # A genotype-linked Number=G field (e.g. PL) whose length disagrees with the genotypes
    # makes `--trim-alt-alleles` abort. If any are present, surgically null just the
    # inconsistent records first (see strip_stale_format; valid likelihoods are kept), then
    # trim and select.
    src = inp
    tmp = None
    try:
        if any(t in format_tags(inp) for t in GENOTYPE_LINKED_FORMAT):
            tmp = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
            strip_stale_format(inp, tmp, fields=GENOTYPE_LINKED_FORMAT, mode="mismatch")
            src = tmp
        sh(f"bcftools view --trim-alt-alleles {q(src)} -Ou "
           f"| bcftools view {select} -{tail}{atomize}", tools=("bcftools",))
    finally:
        if tmp:
            for suffix in ("", ".csi"):
                try:
                    os.remove(tmp + suffix)
                except OSError:
                    pass


class _WhitelistTally:
    """Whether a whitelist has done anything yet, so "it rescued nothing" can be a verdict.

    A whitelist that matches nothing -- wrong contig names, or 1-based positions written into a
    BED -- is worth flagging, since it is otherwise indistinguishable from one that worked. A
    whitelist that simply had nothing to do at one particular step is not: over a chain of
    filters most steps are that, and warning at each of them reads as a misconfiguration.
    """

    def __init__(self) -> None:
        self.deferring = False
        self.rescued = 0
        self.misses: list[str] = []


_WHITELIST = _WhitelistTally()


@contextlib.contextmanager
def deferred_whitelist_warnings():
    """Hold the "rescued nothing" warning until the end of a chain of filters.

    Inside this, a step that rescues nothing is recorded rather than reported, and the warning
    is printed once on the way out -- only if no step in the whole run rescued anything, which
    is the case it was written for.
    """
    prev = (_WHITELIST.deferring, _WHITELIST.rescued, _WHITELIST.misses)
    _WHITELIST.deferring, _WHITELIST.rescued, _WHITELIST.misses = True, 0, []
    try:
        yield _WHITELIST
    finally:
        if not _WHITELIST.rescued and _WHITELIST.misses:
            steps = ", ".join(_WHITELIST.misses)
            print(f"\nWARNING: the whitelist rescued nothing -- no variant it covers would "
                  f"have been dropped by any of the {len(_WHITELIST.misses)} step(s) it "
                  f"applies to ({steps}). Check the contig names, and that its positions are "
                  f"0-based half-open like any BED.")
        _WHITELIST.deferring, _WHITELIST.rescued, _WHITELIST.misses = prev


def _note_whitelist(label: str, keep_bed: str, n: int, dropped_by: str) -> None:
    """Report what the whitelist did at one step: kept variants now, misses at the end."""
    _WHITELIST.rescued += n
    if n:
        say(f"     {label}: {n:,} variant(s) kept by the whitelist that {dropped_by} "
              f"would have dropped")
    elif _WHITELIST.deferring:
        _WHITELIST.misses.append(label)
    else:
        say(f"     {label}: WARNING the whitelist {keep_bed} rescued nothing -- no variant "
              f"it covers would have been dropped here. Check the contig names, and that its "
              f"positions are 0-based half-open like any BED.")


def _rescue_whitelisted(prepared: str, out: str, keep_bed: str | None, label: str) -> int:
    """Add back the whitelisted records a just-run filter dropped. Returns how many.

    The region filters fold a whitelist into their BED, but a filter that judges a record on
    its own numbers -- a QC metric, an allele count, a missingness rate, an allele frequency --
    has no region to carve out. So the whitelist is applied by difference: whatever this step
    dropped that the whitelist covers goes back. The comparison is on position and alleles
    (``bcftools isec``), not on overlap, so a rescued record is that record and not a neighbour.

    ``prepared`` is the input *after* any header-changing part of the step (the FILTER
    annotation, the added FORMAT/ADS, the filled INFO tags) and before its selection. Rescuing
    from the raw input instead would produce records whose header disagrees with the output's,
    which ``bcftools concat`` refuses.

    A rescued record is exempt from this one rule and nothing else: it still faces every later
    step, and it keeps whatever FILTER and INFO the step gave it, so a variant kept despite
    failing QC still says so.
    """
    if not keep_bed:
        return 0
    fmt = out_flag(out)
    tmp = tempfile.mkdtemp()
    wl = os.path.join(tmp, "whitelisted.bcf")
    kept = os.path.join(tmp, "kept.bcf")
    resc = os.path.join(tmp, "rescued.bcf")
    try:
        sh(f"bcftools view {q(prepared)} "
           f"| bedtools intersect -header -a stdin -b {q(keep_bed)} "
           f"| bcftools view -Ob -o {q(wl)}", tools=("bcftools", "bedtools"))
        sh(f"bcftools view {q(out)} -Ob -o {q(kept)}", tools=("bcftools",))
        index_vcf(wl)
        index_vcf(kept)
        sh(f"bcftools isec -C {q(wl)} {q(kept)} -w1 -Ob -o {q(resc)}", tools=("bcftools",))
        n = count_variants(resc)
        if n:
            index_vcf(resc)
            sh(f"bcftools concat -a {q(kept)} {q(resc)} -Ou "
               f"| bcftools sort -O{fmt} -o {q(out)}", tools=("bcftools",))
        _note_whitelist(label, keep_bed, n, "this filter")
        return n
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

def _sorted_bed(bed: str) -> str:
    tmp = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False).name
    sh(f"sort -k1,1 -k2,2n {q(bed)} > {q(tmp)}")
    return tmp


def _whitelisted_record_spans(inp: str, keep_bed: str) -> str:
    """BED of the full span of every input record the whitelist touches, one line per record.

    Record spans, not the whitelist's own intervals, because bedtools sizes a VCF record as
    ``[POS-1, POS-1 + len(REF))``. A multi-base REF therefore reaches past a one-base whitelist
    entry, and carving only that base out of a mask would leave the record still overlapping
    what remained -- it would be dropped despite being whitelisted. Widening to the record's
    own span makes the whitelist mean "keep this variant", which is the useful reading: name a
    position and the variant there survives, whatever its REF length.
    """
    out = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False).name
    # the header has to reach bedtools, which needs it to recognise VCF on stdin
    sh(f"bcftools view {q(inp)} "
       f"| bedtools intersect -a stdin -b {q(keep_bed)} "
       f"| awk -F'\\t' '!/^#/{{print $1\"\\t\"($2-1)\"\\t\"($2-1+length($4))}}' "
       f"| sort -k1,1 -k2,2n > {q(out)}", tools=("bcftools", "bedtools"))
    return out


def _effective_region_bed(inp: str, bed: str, keep_bed: str,
                          exclude: bool) -> tuple[str, int]:
    """Fold a whitelist into a region BED: (effective_bed, n_variants_rescued).

    Applied to the regions rather than by splitting and re-merging the VCF, so the filter stays
    a single bedtools pass:

    * a drop mask has the whitelisted records' spans carved out (``bedtools subtract``), so a
      variant in a whitelisted stretch of a tandem repeat no longer overlaps the mask at all;
    * a keep mask gains them (union), so a variant outside the core survives if whitelisted.

    The count is of records the plain rule would have dropped, which is what tells a caller
    whether the whitelist did anything.
    """
    a, b = _sorted_bed(bed), _whitelisted_record_spans(inp, keep_bed)
    eff = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False).name
    try:
        if exclude:
            sh(f"bedtools subtract -a {q(a)} -b {q(b)} > {q(eff)}", tools=("bedtools",))
            # rescued: whitelisted records that do overlap the mask
            n = _count_bed_lines(f"bedtools intersect -u -a {q(b)} -b {q(a)}")
        else:
            sh(f"cat {q(a)} {q(b)} | sort -k1,1 -k2,2n | bedtools merge > {q(eff)}",
               tools=("bedtools",))
            # rescued: whitelisted records outside the keep regions
            n = _count_bed_lines(f"bedtools intersect -v -a {q(b)} -b {q(a)}")
    finally:
        for f in (a, b):
            if os.path.exists(f):
                os.unlink(f)
    return eff, n


def _count_bed_lines(cmd: str) -> int:
    p = subprocess.run(f"{cmd} | wc -l", shell=True, executable="/bin/bash",
                       stdout=subprocess.PIPE, text=True)
    return int(p.stdout.strip() or 0)


def region_filter(inp: str, out: str, *, bed: str, exclude: bool,
                  keep_bed: str | None = None, label: str = "region_filter") -> int:
    """Keep (``exclude=False``) or drop (``exclude=True``) variants overlapping ``bed``.

    Backs the core-genome (keep), tandem-repeat and paralog (drop) masks. bedtools emits VCF
    text which is piped back through bcftools so any output format works.

    ``keep_bed`` is a whitelist of regions to keep whatever this filter says -- a handful of
    positions inside a tandem repeat that are known to be real, say. A variant overlapping the
    whitelist survives *this* region rule and nothing else: it still faces every other filter
    in the chain.

    Returns the number of variants the whitelist rescued, and says so, since a whitelist that
    matches nothing (wrong contig name, or 1-based positions written into a 0-based BED) is
    otherwise indistinguishable from one that worked. Rescuing nothing *here* is only worth a
    warning if nothing was rescued anywhere -- see :func:`deferred_whitelist_warnings`.
    """
    fmt = out_flag(out)
    v = "-v " if exclude else ""
    eff_bed, rescued = bed, 0
    try:
        if keep_bed:
            eff_bed, rescued = _effective_region_bed(inp, bed, keep_bed, exclude)
            _note_whitelist(label, keep_bed, rescued, "the region rule")
        sh(f"bcftools view {q(inp)} "
           f"| bedtools intersect {v}-header -a stdin -b {q(eff_bed)} "
           f"| bcftools view -O{fmt} -o {q(out)}", tools=("bcftools", "bedtools"))
    finally:
        if keep_bed and os.path.exists(eff_bed):
            os.unlink(eff_bed)
    return rescued


def tandem_repeat_mask(inp: str, out: str, *, bed: str, keep_bed: str | None = None) -> int:
    """Remove variants overlapping a tandem-repeat BED (artifact-prone regions)."""
    return region_filter(inp, out, bed=bed, exclude=True, keep_bed=keep_bed,
                         label="tandem_repeat_mask")


def paralog_mask(inp: str, out: str, *, bed: str, keep_bed: str | None = None) -> int:
    """Remove variants overlapping paralogous / multigene-family genes (mismapping-prone)."""
    return region_filter(inp, out, bed=bed, exclude=True, keep_bed=keep_bed,
                         label="paralog_mask")


def core_region_filter(inp: str, out: str, *, bed: str, keep_bed: str | None = None) -> int:
    """Keep only variants inside the core-genome BED (drop subtelomeric/hypervariable)."""
    return region_filter(inp, out, bed=bed, exclude=False, keep_bed=keep_bed,
                         label="core_region_filter")


def sample_coverage_filter(inp: str, out: str, *, ads_min: int = 10,
                           frac_min: float = 0.80,
                           dropped_samples_path: str | None = None,
                           cov_table_path: str | None = None) -> list[str]:
    """Drop samples covered (ADS >= ads_min) at < frac_min of loci; refresh AC/AN/AF.

    ``cov_table_path`` writes the per-sample coverage table this decision is made from --
    see :func:`sample_coverage_table`. Worth keeping: a dropped sample is otherwise just a
    name in a log, with no way to tell a genuinely thin sample from one that missed the
    threshold by a hair.

    Returns the list of dropped sample names.
    """
    require("bcftools")
    rows = sample_coverage_table(inp, ads_min=ads_min, frac_min=frac_min)
    dropped = sorted(r["sample"] for r in rows if r["dropped"])

    if cov_table_path:
        write_sample_coverage_table(rows, cov_table_path)

    # the borderline cases, said out loud -- an unexpected drop is usually one of these
    for r in rows:
        if r["dropped"] or abs(r["margin"]) <= 0.05:
            detail(f"       {r['sample']}\t{r['n_covered']:,}/{r['n_loci']:,} loci"
                  f"\t{r['frac_covered']:.3f}\tmean ADS {r['mean_ads']:g}"
                  f"\t{'DROPPED' if r['dropped'] else 'kept'}"
                  f" ({r['margin']:+.3f} vs {frac_min:g})")

    if dropped_samples_path:
        with open(dropped_samples_path, "w") as fh:
            fh.write("\n".join(dropped) + ("\n" if dropped else ""))

    fmt = out_flag(out)
    # The singleton re-filter always runs: re-genotyping upstream can leave sites
    # supported by a single sample even when no samples are dropped here.
    keep = 'COUNT(GT!="RR" & GT!="mis") > 1'
    drop_arg = (dropped_samples_path or _write_tmp_list(dropped)) if dropped else None
    view_in = (f"bcftools view -S ^{q(drop_arg)} {q(inp)} -Ou"
               if drop_arg else f"bcftools view {q(inp)} -Ou")
    sh(f"{view_in} "
       f"| bcftools view -i {q(keep)} -Ou "
       f"| bcftools +fill-tags -O{fmt} -o {q(out)} -- -t AC,AN,AF", tools=("bcftools",))
    return dropped


def locus_missingness_filter(inp: str, out: str, *, f_missing_max: float = 0.05,
                             ads_min: int = 10, sample_frac_min: float = 0.95,
                             keep_bed: str | None = None) -> int:
    """Keep loci with < f_missing_max missing AND >= sample_frac_min at ADS >= ads_min.

    ``keep_bed`` whitelists regions from this rule, for a locus worth keeping even where it
    is thinly covered.
    """
    fmt = out_flag(out)
    expr = (f"F_MISSING < {f_missing_max} & "
            f"COUNT(FMT/ADS>={ads_min})/N_SAMPLES >= {sample_frac_min}")
    tags = f"bcftools annotate -x INFO/F_MISSING {q(inp)} -Ou | bcftools +fill-tags"
    if not keep_bed:
        sh(f"{tags} -Ou -- -t AC,AN,AF,F_MISSING "
           f"| bcftools view -i {q(expr)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return 0
    prep = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
    try:
        sh(f"{tags} -Ob -o {q(prep)} -- -t AC,AN,AF,F_MISSING", tools=("bcftools",))
        sh(f"bcftools view -i {q(expr)} {q(prep)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return _rescue_whitelisted(prep, out, keep_bed, "locus_missingness_filter")
    finally:
        if os.path.exists(prep):
            os.unlink(prep)


def _maf_expr(maf_min: float, maf_max: float | None) -> str:
    """The frequency window as a bcftools expression on ``INFO/MAF``.

    ``MAF`` (from ``+fill-tags``) is the frequency of the **second most common** allele --
    the textbook minor-allele frequency, and the same quantity the grouped path computes, so
    the two cannot mean different things.

    Neither of the obvious ``bcftools view`` shortcuts says this:

    * ``-q`` defaults to ``nref``, the **sum** of every alternate's frequency. Three
      alternates at 1% each sum past a 2% floor, so a site whose every allele is below the
      floor survives it.
    * ``-q X:minor`` is the **least frequent** allele, not the second most common. On a site
      at 0.60 / 0.35 / 0.05 that is 0.05, so a 0.10 floor drops a site with a 35% minor
      allele. It is a far stricter rule than the name suggests, and on a multiallelic
      callset it removes almost everything.

    On a biallelic site all three coincide, which is why this only surfaced once multiallelic
    SNPs could reach the step.

    The two bounds ask different questions, so they are spelled differently:

    * ``maf_max`` unset (the usual case) is a **minor-allele floor**: ``MAF >= maf_min``. A
      floor on the minor allele already implies a ``1 - maf_min`` ceiling on the major one,
      so the symmetric window needs nothing else.
    * an explicit ``maf_max`` is an **alternate-frequency band**, which is what asking for an
      asymmetric window means: ``maf_min <= MAX(AF) <= maf_max``, on the most common
      alternate. Reading the floor as a minor-allele floor here would quietly make an
      asymmetric request stricter than it looks -- a site at AF 0.75 has a 0.25 minor allele,
      and would fail a 0.3 floor despite sitting inside a 0.3-0.8 band.

    Both spellings select exactly what the old ``-q``/``-Q`` pair did on a biallelic site.
    """
    if maf_max is None:
        return f"MAF >= {maf_min}"
    return f"MAX(AF) >= {maf_min} && MAX(AF) <= {maf_max}"


def maf_filter(inp: str, out: str, *, maf_min: float = 0.01, maf_max: float | None = None,
               meta: str | None = None, group_col: str | None = None,
               sample_col: str = "sample", keep_bed: str | None = None) -> int:
    """Drop rare and near-fixed alleles by an allele-frequency window ``[maf_min, maf_max]``.

    The two bounds are usually symmetric (a 0.02 floor pairs with a 0.98 ceiling), so
    ``maf_max`` defaults to ``1 - maf_min`` when left unset; pass it explicitly for an
    asymmetric window. Left unset the test is a minor-allele floor (``INFO/MAF``, the
    frequency of the second most common allele); given a value it becomes a band on the most
    common alternate, which is what an asymmetric request means. Either way the window reads
    the same however many alleles a site has -- see :func:`_maf_expr`.

    With ``meta`` + ``group_col`` (a per-sample metadata table with ``sample_col`` and
    ``group_col`` columns), the frequency is judged **per group** and a site is kept if its
    minor-allele frequency (``+fill-tags -t MAF``, the same quantity the ungrouped
    window tests) is >= ``maf_min`` in **any** group — computed on the combined VCF
    (never split-and-merged), so every sample's genotypes are preserved even at a site that
    is monomorphic in its own group but polymorphic elsewhere. ``maf_max`` is not used in
    grouped mode (the criterion is a per-group minor-allele-frequency floor).

    ``keep_bed`` whitelists regions from the frequency window. This is the one most worth
    reaching for: a resistance allele can sit at a few percent in one cohort and still be the
    thing being looked for, and a MAF floor is exactly what removes it.
    """
    if meta and group_col:
        return _maf_filter_grouped(inp, out, meta=meta, group_col=group_col,
                                   sample_col=sample_col, maf_min=maf_min, keep_bed=keep_bed)
    fmt = out_flag(out)
    expr = _maf_expr(maf_min, maf_max)
    if not keep_bed:
        sh(f"bcftools +fill-tags {q(inp)} -Ou -- -t AC,AN,AF,MAF "
           f"| bcftools view -i {q(expr)} -O{fmt} -o {q(out)}",
           tools=("bcftools",))
        return 0
    prep = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
    try:
        sh(f"bcftools +fill-tags {q(inp)} -Ob -o {q(prep)} -- -t AC,AN,AF,MAF",
           tools=("bcftools",))
        sh(f"bcftools view -i {q(expr)} {q(prep)} -O{fmt} -o {q(out)}",
           tools=("bcftools",))
        return _rescue_whitelisted(prep, out, keep_bed, "maf_filter")
    finally:
        if os.path.exists(prep):
            os.unlink(prep)


def _vcf_samples(path: str) -> set[str]:
    p = subprocess.run(["bcftools", "query", "-l", str(path)],
                       stdout=subprocess.PIPE, text=True)
    return {s for s in p.stdout.splitlines() if s}


def _maf_filter_grouped(inp: str, out: str, *, meta: str, group_col: str,
                        sample_col: str, maf_min: float,
                        keep_bed: str | None = None) -> int:
    require("bcftools")
    present = _vcf_samples(inp)
    groups: dict[str, list[str]] = defaultdict(list)
    with open(meta) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        # `Sample` and `sample` are the same column to everyone but a string comparison
        s_col = Utils.resolve_column(fields, sample_col, source=f"metadata ({meta})")
        g_col = Utils.resolve_column(fields, group_col, source=f"metadata ({meta})")
        for got, want in ((s_col, sample_col), (g_col, group_col)):
            if got != want:
                print(f"  note: metadata column '{got}' read as '{want}'")
        for row in reader:
            s, g = row[s_col], row[g_col]
            if s in present and g:
                groups[g].append(s)
    if not groups:
        raise SystemExit(f"ERROR: no samples in {meta} overlap the VCF")

    tmp = tempfile.mkdtemp(prefix="maf_grouped_")
    try:
        union = os.path.join(tmp, "pass_positions.txt")
        with open(union, "w") as acc:
            for g, samples in groups.items():
                sfile = os.path.join(tmp, "samples.txt")
                with open(sfile, "w") as sf:
                    sf.write("\n".join(samples) + "\n")
                # per-group minor-allele frequency, keep sites with MAF >= maf_min
                p = subprocess.run(
                    f"bcftools view -S {q(sfile)} --force-samples {q(inp)} -Ou "
                    f"| bcftools +fill-tags -Ou -- -t MAF "
                    f"| bcftools query -i 'MAF>={maf_min}' -f '%CHROM\\t%POS\\n'",
                    shell=True, executable="/bin/bash",
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if p.returncode != 0:
                    raise SystemExit(f"ERROR: per-group MAF failed for '{g}':\n{p.stderr}")
                acc.write(p.stdout)
        sorted_pos = os.path.join(tmp, "pass_sorted.txt")
        sh(f"LC_ALL=C sort -k1,1 -k2,2n -u {q(union)} > {q(sorted_pos)}")
        fmt = out_flag(out)
        # apply the union of passing sites to the ORIGINAL combined VCF (all samples kept)
        sh(f"bcftools view -T {q(sorted_pos)} {q(inp)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        # selection is by position list off the untouched input, so its header already
        # matches and the whitelist can be rescued straight from it
        rescued = _rescue_whitelisted(inp, out, keep_bed, "maf_filter")
    finally:
        for f in os.listdir(tmp):
            try:
                os.remove(os.path.join(tmp, f))
            except OSError:
                pass
        try:
            os.rmdir(tmp)
        except OSError:
            pass
    return rescued


def vcf_to_bed(inp: str, out: str | None = None, *, snps_only: bool = False,
               name_column: bool = True) -> None:
    """Write the records of a VCF/BCF as BED. ``out=None`` writes to stdout.

    Columns are ``chrom``, 0-based start, end, and (unless ``name_column=False``) a
    ``chrom:pos0`` name — the canonical label
    (:func:`~plasgenomicsutils.lib.intervals.snp_label`) a panel loaded from either BED or
    VCF derives, so the two agree.

    **Everything is 0-based half-open**, which is the whole point of the conversion: VCF
    ``POS`` is 1-based and BED is not, and doing this by hand is where the off-by-one
    lives. The interval is the span of the REF allele — one base for a SNP, ``len(REF)``
    for an indel or MNP — so a record's extent is what a region file needs it to be rather
    than just its start.

    ``snps_only`` keeps single-base substitutions and drops indels and everything else.
    """
    require("bcftools")
    fields = "%CHROM\\t%POS0\\t%END" + ("\\t%CHROM:%POS0" if name_column else "") + "\\n"
    pipe = f"bcftools view -v snps {q(inp)} -Ou | " if snps_only else ""
    src = "-" if snps_only else q(inp)
    redirect = f" > {q(out)}" if out else ""
    sh(f"{pipe}bcftools query -f '{fields}' {src}{redirect}", tools=("bcftools",))


def snp_bed(inp: str, bed: str) -> None:
    """Write a BED of the SNP positions in ``inp`` — the SNP panel the IBD tools read
    (``build_ibd_matrix --snp-format bed``).

    A :func:`vcf_to_bed` restricted to SNPs; see there for the columns and the coordinate
    convention.
    """
    vcf_to_bed(inp, bed, snps_only=True)


# --- small bcftools query helpers -------------------------------------------

def _count(path: str) -> int:
    p = subprocess.run(f"bcftools view -H {q(path)} | wc -l", shell=True,
                       executable="/bin/bash", stdout=subprocess.PIPE, text=True)
    return int(p.stdout.strip() or 0)


def _samples(path: str) -> list[str]:
    p = subprocess.run(f"bcftools query -l {q(path)}", shell=True,
                       executable="/bin/bash", stdout=subprocess.PIPE, text=True)
    return [s for s in p.stdout.splitlines() if s]


def _has_format_tag(path: str, tag: str) -> bool:
    p = subprocess.run(f"bcftools view -h {q(path)}", shell=True, executable="/bin/bash",
                       stdout=subprocess.PIPE, text=True)
    return f"##FORMAT=<ID={tag}," in p.stdout


def _per_sample_coverage(path: str, ads_min: int) -> dict[str, dict[str, int]]:
    """Per-sample ADS summary in one pass: loci seen, loci covered, missing, total depth.

    ADS is a scalar (``int(smpl_sum(FORMAT/AD))``), so testing ``>= ads_min`` here selects
    exactly what ``bcftools view -i 'FMT/ADS>=n'`` selects. That matters: the coverage table
    and the keep/drop decision are both read off this one result, so the table always accounts
    for the decision rather than being a second measurement that could disagree with it.

    Aggregated in awk, since a cohort-sized callset is millions of sample-by-site lines.
    """
    prog = (
        'BEGIN{FS="\\t"} '
        '{n[$1]++; if($2=="."||$2==""){miss[$1]++} '
        'else{s[$1]+=$2; if($2+0>=MIN){cov[$1]++}}} '
        'END{for(k in n) printf "%s\\t%d\\t%d\\t%d\\t%d\\n", k, n[k], cov[k]+0, miss[k]+0,'
        ' s[k]+0}'
    )
    cmd = (f"bcftools query -f '[%SAMPLE\\t%ADS\\n]' {q(path)} "
           f"| awk -v MIN={int(ads_min)} {q(prog)}")
    p = subprocess.run(cmd, shell=True, executable="/bin/bash",
                       stdout=subprocess.PIPE, text=True)
    out: dict[str, dict[str, int]] = {}
    for line in p.stdout.splitlines():
        if not line.strip():
            continue
        # sample names can contain tabs in no sane file, but split from the right anyway
        name, n_seen, n_cov, n_miss, ads_sum = line.rsplit("\t", 4)
        out[name] = {"n_loci": int(n_seen), "n_covered": int(n_cov),
                     "n_missing_ads": int(n_miss), "ads_sum": int(ads_sum)}
    return out


def sample_coverage_table(path: str, *, ads_min: int = 10,
                          frac_min: float = 0.80) -> list[dict]:
    """Per-sample coverage, and whether :func:`sample_coverage_filter` drops each sample.

    One row per sample, worst-covered first, so a surprising drop is the first thing read.
    ``margin`` is ``frac_covered - frac_min``: negative means dropped, and a value near zero
    means the sample sat on the threshold rather than being obviously bad. ``mean_ads`` and
    ``n_missing_ads`` say *why* -- thin coverage everywhere reads differently from a sample
    that is simply absent at many sites.
    """
    # Without FORMAT/ADS every sample scores zero covered loci and the filter quietly drops
    # the entire cohort. bcftools does say so on stderr, but the counts come back empty rather
    # than failing, so refuse instead -- silently emptying a callset is the worst outcome here.
    if not _has_format_tag(path, "ADS"):
        raise SystemExit(
            f"ERROR: {path} has no FORMAT/ADS, which is what sample coverage is measured on. "
            "Run singleton_filter_add_ads first (it adds ADS as int(smpl_sum(FORMAT/AD))), or "
            "add the tag with: bcftools +fill-tags -- -t "
            "'FORMAT/ADS=int(smpl_sum(FORMAT/AD))'.")
    stats = _per_sample_coverage(path, ads_min)
    total = _count(path)
    rows = []
    for name in _samples(path):
        st = stats.get(name, {"n_loci": 0, "n_covered": 0, "n_missing_ads": 0, "ads_sum": 0})
        n_loci = st["n_loci"] or total
        frac = (st["n_covered"] / n_loci) if n_loci else 0.0
        with_ads = n_loci - st["n_missing_ads"]
        rows.append({
            "sample": name,
            "n_loci": n_loci,
            "n_covered": st["n_covered"],
            "frac_covered": round(frac, 6),
            "n_missing_ads": st["n_missing_ads"],
            "mean_ads": round(st["ads_sum"] / with_ads, 2) if with_ads else 0.0,
            "ads_min": ads_min,
            "frac_min": frac_min,
            "margin": round(frac - frac_min, 6),
            "dropped": frac < frac_min,
        })
    return sorted(rows, key=lambda r: r["frac_covered"])


COVERAGE_TABLE_COLUMNS = ["sample", "n_loci", "n_covered", "frac_covered", "n_missing_ads",
                          "mean_ads", "ads_min", "frac_min", "margin", "dropped"]


def write_sample_coverage_table(rows: list[dict], path: str) -> None:
    """Write the per-sample coverage table that explains each keep/drop decision."""
    with open(path, "w") as fh:
        fh.write("\t".join(COVERAGE_TABLE_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COVERAGE_TABLE_COLUMNS) + "\n")


def _write_tmp_list(names: list[str]) -> str:
    import tempfile
    fh = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    fh.write("\n".join(names) + "\n")
    fh.close()
    return fh.name
