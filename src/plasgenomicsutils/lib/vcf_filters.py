"""Parameterized VCF/BCF filtering steps, each backed by bcftools/bedtools.

Every step takes an input and output path plus tunable thresholds, so the same
building blocks can be run individually or chained by ``filter_pipeline``.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

from ..utils.small_utils import Utils
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


def _bcftools_qc_expr(*, qd, mq, strand_bias_p, max_bias_z, read_pos_z, mq0f, bqbz_z):
    """The bcftools-native mirror of the GATK hard filter, as a bcftools -e expression.

    Each GATK metric has a counterpart in what `bcftools mpileup` writes, but two of them
    do not carry over naively:

    * GATK ``FS`` is Phred-scaled and rises with bias; **bcftools ``FS`` is the p-value
      itself**, so the test is ``FS <`` a small number rather than ``FS >`` a large one.
      (Phred 60, GATK's usual cutoff, is p = 1e-6 -- the same statement either way.)
    * GATK's rank-sum filters are one-sided, because the sign says which way the alt
      allele leans. The bcftools ``*BZ`` tags are documented as "closer to 0 is better",
      and a read-position artifact shows up as either sign, so they are tested on
      ``abs()``.

    ``DP`` is written as ``INFO/DP`` deliberately: a bare ``DP`` in a bcftools expression
    resolves to ``FORMAT/DP`` where both exist, which silently changes what is filtered.
    """
    parts = []
    if qd is not None:
        parts.append(f"QUAL/INFO/DP < {qd}")
    if mq is not None:
        parts.append(f"MQ < {mq}")
    if strand_bias_p is not None:                       # strand bias, Fisher p-value
        parts.append(f'(FS!="." && FS < {strand_bias_p})')
    if read_pos_z is not None:                          # the alt sits near the read ends
        parts.append(f'(RPBZ!="." && abs(RPBZ) > {read_pos_z})')
        parts.append(f'(SCBZ!="." && abs(SCBZ) > {read_pos_z})')
    if max_bias_z is not None:                          # mapping quality, and vs strand
        parts.append(f'(MQBZ!="." && abs(MQBZ) > {max_bias_z})')
        parts.append(f'(MQSBZ!="." && abs(MQSBZ) > {max_bias_z})')
    if bqbz_z is not None:
        parts.append(f'(BQBZ!="." && abs(BQBZ) > {bqbz_z})')
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
    "max_bias_z": ("MQBZ", "MQSBZ"), "bqbz_z": ("BQBZ",), "mq0f": ("MQ0F",),
    "mq": ("MQ",),
}


def _check_bcftools_qc_tags(inp: str, wanted: dict) -> None:
    have = info_tags(inp)
    missing = sorted({t for k, v in wanted.items() if v is not None
                      for t in _BCFTOOLS_QC_NEEDS[k] if t not in have})
    if not missing:
        return
    raise SystemExit(
        "hard_qc_filter --caller bcftools: this callset has no INFO/"
        + ", INFO/".join(missing) + ".\n"
        "  A comparison against a tag that is not there is simply false, so the filter "
        "would keep\n  everything and say nothing. Call with the annotations it needs:\n\n"
        f"    {BCFTOOLS_CALL_RECIPE}\n\n"
        "  ...or disable the thresholds that read the missing tags."
    )


def no_alt_filter(inp: str, out: str, *, keep: bool = False,
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
    """
    n_no_alt = count_variants_matching(inp, 'ALT="."')
    fmt = out_flag(out)
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
                   mq: float | None = 55, sor: float = 3, mqranksum: float = -5.0,
                   readposranksum: float = -5.0, fs: float | None = None,
                   strand_bias_p: float | None = 1e-6, read_pos_z: float | None = 5.0,
                   max_bias_z: float | None = 5.0, bqbz_z: float | None = None,
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

    if caller == "bcftools":
        thresholds = dict(qd=qd, mq=mq, strand_bias_p=strand_bias_p,
                          read_pos_z=read_pos_z, max_bias_z=max_bias_z, bqbz_z=bqbz_z,
                          mq0f=mq0f)
        _check_bcftools_qc_tags(inp, thresholds)
        expr = _bcftools_qc_expr(**thresholds)
    else:
        parts = [
            f"QD < {qd}",
            f"MQ < {mq}",
            f"SOR > {sor}",
            f'(MQRankSum!="." && MQRankSum < {mqranksum})',
            f'(ReadPosRankSum!="." && ReadPosRankSum < {readposranksum})',
        ]
        if fs is not None:
            parts.append(f"FS > {fs}")
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


def biallelic_snp_filter(inp: str, out: str, *, trim: bool = True) -> None:
    """Keep biallelic SNPs, optionally trimming unused ALT alleles first.

    When ``trim`` is set (default), ALT alleles absent from every genotype are
    dropped before the biallelic test (``bcftools view --trim-alt-alleles``). Run
    this AFTER re-genotyping (``filter_ad_regenotype``): a site that looked
    multiallelic only because one sample carried a low-level artifact allele —
    now re-genotyped away — collapses to a genuine biallelic SNP and is kept,
    instead of being discarded by a naive ``-m2 -M2``. Sites left with no ALT
    after trimming are dropped by ``-v snps``.
    """
    fmt = out_flag(out)
    if not trim:
        sh(f"bcftools view -m2 -M2 -v snps {q(inp)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return
    # A genotype-linked Number=G field (e.g. PL) whose length disagrees with the genotypes
    # makes `--trim-alt-alleles` abort. If any are present, surgically null just the
    # inconsistent records first (see strip_stale_format; valid likelihoods are kept), then
    # trim and keep biallelic SNPs.
    if any(t in format_tags(inp) for t in GENOTYPE_LINKED_FORMAT):
        tmp = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
        try:
            strip_stale_format(inp, tmp, fields=GENOTYPE_LINKED_FORMAT, mode="mismatch")
            sh(f"bcftools view --trim-alt-alleles {q(tmp)} -Ou "
               f"| bcftools view -m2 -M2 -v snps - -O{fmt} -o {q(out)}", tools=("bcftools",))
        finally:
            for suffix in ("", ".csi"):
                try:
                    os.remove(tmp + suffix)
                except OSError:
                    pass
    else:
        sh(f"bcftools view --trim-alt-alleles {q(inp)} -Ou "
           f"| bcftools view -m2 -M2 -v snps - -O{fmt} -o {q(out)}", tools=("bcftools",))



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
            print(f"     {label}: {n:,} variant(s) kept by the whitelist that this filter "
                  f"would have dropped")
        else:
            print(f"     {label}: WARNING the whitelist {keep_bed} rescued nothing -- no "
                  f"variant it covers would have been dropped here. Check the contig names, "
                  f"and that its positions are 0-based half-open like any BED.")
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
    otherwise indistinguishable from one that worked.
    """
    fmt = out_flag(out)
    v = "-v " if exclude else ""
    eff_bed, rescued = bed, 0
    try:
        if keep_bed:
            eff_bed, rescued = _effective_region_bed(inp, bed, keep_bed, exclude)
            if rescued:
                print(f"     {label}: {rescued:,} variant(s) kept by the whitelist that the "
                      f"region rule would have dropped")
            else:
                print(f"     {label}: WARNING the whitelist {keep_bed} rescued nothing -- "
                      f"no variant it covers would have been dropped here. Check the contig "
                      f"names, and that its positions are 0-based half-open like any BED.")
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
            print(f"       {r['sample']}\t{r['n_covered']:,}/{r['n_loci']:,} loci"
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


def maf_filter(inp: str, out: str, *, maf_min: float = 0.01, maf_max: float | None = None,
               meta: str | None = None, group_col: str | None = None,
               sample_col: str = "sample", keep_bed: str | None = None) -> int:
    """Drop rare and near-fixed alleles by an allele-frequency window ``[maf_min, maf_max]``.

    The two bounds are usually symmetric (a 0.02 floor pairs with a 0.98 ceiling), so
    ``maf_max`` defaults to ``1 - maf_min`` when left unset; pass it explicitly for an
    asymmetric window.

    With ``meta`` + ``group_col`` (a per-sample metadata table with ``sample_col`` and
    ``group_col`` columns), the frequency is judged **per group** and a site is kept if its
    minor-allele frequency is >= ``maf_min`` in **any** group — computed on the combined VCF
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
    if maf_max is None:
        maf_max = 1 - maf_min
    fmt = out_flag(out)
    if not keep_bed:
        sh(f"bcftools +fill-tags {q(inp)} -Ou -- -t AC,AN,AF "
           f"| bcftools view -q {maf_min} -Q {maf_max} -O{fmt} -o {q(out)}",
           tools=("bcftools",))
        return 0
    prep = tempfile.NamedTemporaryFile(suffix=".bcf", delete=False).name
    try:
        sh(f"bcftools +fill-tags {q(inp)} -Ob -o {q(prep)} -- -t AC,AN,AF", tools=("bcftools",))
        sh(f"bcftools view -q {maf_min} -Q {maf_max} {q(prep)} -O{fmt} -o {q(out)}",
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


def snp_bed(inp: str, bed: str) -> None:
    """Write a BED of the SNP positions in ``inp`` — the SNP panel the IBD tools read
    (``build_ibd_matrix --snp-format bed``).

    Columns: ``chrom``, 0-based start, end, and a ``chrom:pos0`` name. Every column is
    0-based, so the file does not contradict itself, and the name matches the canonical
    label (:func:`~plasgenomicsutils.lib.intervals.snp_label`) that a panel loaded from
    either BED or VCF derives. Only SNP records are emitted.
    """
    require("bcftools")
    sh(f"bcftools view -v snps {q(inp)} -Ou "
       f"| bcftools query -f '%CHROM\\t%POS0\\t%END\\t%CHROM:%POS0\\n' - > {q(bed)}",
       tools=("bcftools",))


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
