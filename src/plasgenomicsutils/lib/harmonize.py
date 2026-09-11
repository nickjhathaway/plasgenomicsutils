"""Harmonize separately-called cohorts so they can be merged.

When variant calling is done separately across cohorts at the same sites, the
per-file ALT sets and allele-count INFO fields disagree, and ``bcftools merge``
then chokes on the mismatched cardinalities. This:

  1. cleans spurious low-level ALTs per file and re-genotypes,
  2. builds the union of real ALTs at each site across all files, and
  3. rewrites each file to that union (padding AD with zeros, re-genotyping only
     samples that gained a new allele),

stripping the stale allele-count INFO fields so ``bcftools +fill-tags`` can
recompute them after the merge. Two streaming passes are used (accumulate the
union, then rewrite), so records are never all held in memory at once. Inputs
must be coordinate-sorted.

Every ``Number=R`` count field moves with the alleles -- FORMAT and INFO AD, ADF, ADR
are re-laid-out on the union with zeros for alleles a file did not carry -- and the
``Number=A``/``G`` INFO fields (AC, AF) are dropped. :mod:`.call_variants` uses the same
two passes with no cleaning and ``regenotype=False`` to join a cohort it called in groups
of alignments; there the files are one caller's output split by sample, and the
genotypes are re-indexed rather than re-called.
"""

from __future__ import annotations

import numpy as np

from .ad_genotype import regenotype_from_ad


# INFO fields whose values depend on the allele set / count and become
# inconsistent once alleles are edited; recomputed downstream by +fill-tags.
STALE_INFO_FIELDS = [
    "AC", "AF", "AN", "DP4", "VDB", "SGB",
    "RPBZ", "MQBZ", "MQSBZ", "BQBZ", "SCBZ", "MQ0F",
]


def strip_stale_info(rec, fields=STALE_INFO_FIELDS) -> None:
    """Drop INFO fields whose value is tied to the allele set.

    ``Number=A/R`` fields (AC, AF, DP4 and friends) are sized by how many alleles the
    record has, so once the ALT set changes they are both wrong and a merge hazard:
    ``bcftools merge`` aborts on an allele-count mismatch. ``MQ`` is kept -- it is a single
    value that does not depend on the alleles. ``fields`` is the list to drop; the default
    is what the cross-cohort harmonize strips, and a caller that reshapes the per-allele
    counts itself (see :func:`field_shapes`) can pass a shorter one.
    """
    for field in fields:
        try:
            del rec.info[field]
        except KeyError:
            pass


def field_shapes(header) -> dict:
    """Which header fields are shaped by the allele list, so they can be reshaped with it.

    Returns ``{"info_r": [...], "format_r": [...], "info_ag": [...]}``: the ``Number=R``
    INFO and FORMAT fields (one value per allele, REF first -- AD, ADF, ADR), which are
    padded with zeros or cut down alongside the alleles, and the ``Number=A``/``G`` INFO
    fields (AC, AF), which cannot be and are dropped. FORMAT ``Number=A``/``G`` fields (PL)
    are left to :func:`stale_format_fields`, since pysam cannot resize them cleanly and
    `bcftools annotate -x` does it afterwards.
    """
    info_r, info_ag, format_r = [], [], []
    for fid, meta in header.info.items():
        if str(meta.number) == "R":
            info_r.append(fid)
        elif str(meta.number) in ("A", "G"):
            info_ag.append(fid)
    for fid, meta in header.formats.items():
        if str(meta.number) == "R":
            format_r.append(fid)
    return {"info_r": info_r, "format_r": format_r, "info_ag": info_ag}


def _snapshot_r(rec, shapes: dict):
    """Copy every ``Number=R`` value off a record before its alleles are changed.

    pysam reinterprets the per-allele FORMAT/INFO arrays the moment ``rec.alleles`` is
    assigned, so the old values have to be read first. A sample with a missing or
    wrong-length array snapshots as ``None``.
    """
    n = len(rec.alleles)
    info = {}
    for fid in shapes["info_r"]:
        v = rec.info.get(fid, None)
        info[fid] = None if (v is None or any(x is None for x in v) or len(v) != n) else list(v)
    fmt = {}
    for sname in rec.samples:
        per = {}
        for fid in shapes["format_r"]:
            v = rec.samples[sname].get(fid, None)
            per[fid] = None if (v is None or any(x is None for x in v) or len(v) != n) \
                else list(v)
        fmt[sname] = per
    return info, fmt


def _remap_gt(sample, old_to_new: dict) -> int:
    """Re-index a genotype through an old-allele -> new-allele map, keeping phasing.

    An allele with no entry in the map is one the union does not carry, which should be
    impossible: the union is built from these same records. It became reachable through the
    pass-1/pass-2 disagreement that :func:`harmonize_file` now refuses outright. Should any
    other route to it appear, the call becomes missing and is **counted** rather than
    raising a bare ``KeyError`` that names neither the position nor the file. Returns how
    many alleles were dropped, so the caller can say so.
    """
    gt = sample.get("GT", None)
    if gt is None:
        return 0
    phased = sample.phased
    new = tuple(None if g is None else old_to_new.get(g) for g in gt)
    dropped = sum(1 for g, n in zip(gt, new) if g is not None and n is None)
    sample["GT"] = new
    sample.phased = phased
    return dropped


def clean_record(rec, min_ad: int, min_af: float, het_min_af: float, *,
                 regenotype: bool = True, shapes: dict | None = None,
                 stale_info=STALE_INFO_FIELDS) -> None:
    """Zero low-level ALT depths, drop empty ALTs, re-genotype (mutates ``rec``).

    Records reduced to ref-only are kept with ALT="." so the union step can
    still fill ALTs from other files. REF depth is never zeroed.

    ``regenotype=False`` keeps the caller's genotypes where it can: when an allele is
    dropped, a sample whose genotype does not use it is re-indexed rather than re-called
    from AD, and only a sample that was carrying the dropped allele is re-called. With
    ``min_ad=0`` and ``min_af=0`` (no cleaning) that leaves every genotype as called and
    only removes alleles no sample has a read for. ``shapes`` (:func:`field_shapes`) lists
    the other per-allele fields to cut down with the alleles; it is read off the header
    when not given. ``stale_info`` is what :func:`strip_stale_info` drops, plus every
    ``Number=A``/``G`` INFO field.
    """
    alleles = list(rec.alleles)
    n_alleles = len(alleles)
    if shapes is None:
        shapes = field_shapes(rec.header)
    stale = list(stale_info) + [f for f in shapes["info_ag"] if f not in stale_info]
    if n_alleles < 2 or (n_alleles == 2 and alleles[1] == "."):
        return  # already ref-only

    total_ad_per_allele = [0] * n_alleles
    new_ad_per_sample = {}
    for sname in rec.samples:
        ad = rec.samples[sname].get("AD", None)
        if ad is None or any(v is None for v in ad):
            new_ad_per_sample[sname] = None
            continue
        ad = list(ad)
        if len(ad) != n_alleles:
            new_ad_per_sample[sname] = None
            continue
        total_depth = sum(ad)
        cleaned = []
        for i, count in enumerate(ad):
            if i == 0:
                cleaned.append(count)  # never zero REF
            else:
                af = count / total_depth if total_depth > 0 else 0.0
                cleaned.append(0 if (count < min_ad or af < min_af) else count)
        new_ad_per_sample[sname] = cleaned
        for i, v in enumerate(cleaned):
            total_ad_per_allele[i] += v

    alts_to_keep = [True] * n_alleles
    any_removed = False
    for i in range(1, n_alleles):
        if total_ad_per_allele[i] == 0:
            alts_to_keep[i] = False
            any_removed = True

    if any_removed:
        keep_idx = [i for i in range(n_alleles) if alts_to_keep[i]]
        old_to_new = {old: new for new, old in enumerate(keep_idx)}
        info_r, fmt_r = _snapshot_r(rec, shapes)
        gts = {sname: rec.samples[sname].get("GT", None) for sname in rec.samples}
        new_alleles = [alleles[i] for i in keep_idx]
        if len(new_alleles) == 1:
            rec.alleles = (alleles[0], ".")
            for sname in rec.samples:
                ad = new_ad_per_sample.get(sname)
                ref_ad = ad[0] if ad is not None else 0
                rec.samples[sname]["AD"] = (ref_ad, 0)
                # Reference only if there are reference reads to say so. This used to stamp
                # 0/0 on every sample, including one with no AD and one with AD=0,0 -- a
                # sample with no reads at all called confidently homozygous reference, which
                # is the rule the rest of the module refuses ("total 0 -> missing") and the
                # gVCF trap the singleton counter documents.
                rec.samples[sname]["GT"] = (0, 0) if ref_ad > 0 else (None, None)
                # the other per-allele counts shrink to REF plus the empty ALT slot
                for fid in shapes["format_r"]:
                    if fid != "AD":
                        v = fmt_r[sname][fid]
                        rec.samples[sname][fid] = (v[0] if v else 0, 0)
            for fid in shapes["info_r"]:
                if info_r[fid] is not None:
                    rec.info[fid] = (info_r[fid][0], 0)
            strip_stale_info(rec, stale)
            return
        n_new = len(new_alleles)
        rec.alleles = tuple(new_alleles)
        for fid in shapes["info_r"]:
            if info_r[fid] is not None:
                rec.info[fid] = tuple(info_r[fid][i] for i in keep_idx)
        for sname in rec.samples:
            sample = rec.samples[sname]
            ad = new_ad_per_sample.get(sname)
            for fid in shapes["format_r"]:
                if fid == "AD":
                    continue
                v = fmt_r[sname][fid]
                sample[fid] = tuple(v[i] for i in keep_idx) if v else tuple([0] * n_new)
            if ad is not None and len(ad) == n_alleles:
                new_ad = [ad[i] for i in keep_idx]
                sample["AD"] = tuple(new_ad)
                gt = gts[sname]
                uses_dropped = gt is not None and any(g is not None and g not in old_to_new
                                                      for g in gt)
                if regenotype or uses_dropped or gt is None:
                    gt = regenotype_from_ad(new_ad, het_min_af)
                    sample["GT"] = gt if gt is not None else (None, None)
                else:
                    _remap_gt(sample, old_to_new)
            else:
                # Sample had no usable AD; still write a correct-length AD so the
                # record stays consistent (Number=R) after the allele count drops.
                sample["AD"] = tuple([0] * n_new)
                sample["GT"] = (None, None)
        strip_stale_info(rec, stale)
    else:
        for sname in rec.samples:
            sample = rec.samples[sname]
            ad = new_ad_per_sample.get(sname)
            if ad is None or len(ad) != n_alleles:
                continue
            original_ad = list(sample.get("AD", []) or [])
            sample["AD"] = tuple(ad)
            if ad != original_ad and regenotype:
                gt = regenotype_from_ad(ad, het_min_af)
                sample["GT"] = gt if gt is not None else (None, None)
        strip_stale_info(rec, stale)


def n_real_alts(rec) -> int:
    """Number of non-``.`` ALT alleles on a record."""
    return sum(1 for a in rec.alleles[1:] if a != ".")


def surviving_alleles(rec, min_ad: int, min_af: float):
    """(ref, set-of-surviving-ALTs) after AD cleaning, without mutating ``rec``.

    Pass 1 (union building) only needs which ALTs still have support once
    sub-threshold per-sample depths are zeroed — not the re-genotyped record. This
    computes that with one vectorized numpy pass and no genotype work, returning
    exactly the ALT set :func:`clean_record` would leave on the record.
    """
    alleles = list(rec.alleles)
    n = len(alleles)
    if n < 2 or (n == 2 and alleles[1] == "."):
        return alleles[0], set()

    rows = []
    for sname in rec.samples:
        ad = rec.samples[sname].get("AD", None)
        if ad is None or len(ad) != n or any(v is None for v in ad):
            continue
        rows.append(ad)
    if not rows:
        return alleles[0], set()

    A = np.asarray(rows, dtype=float)              # (m_valid, n_alleles)
    depth = A.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        af = A / depth
    fail = (A < min_ad) | (af < min_af)            # per-sample sub-threshold ALTs
    fail[:, 0] = False                             # REF is never zeroed
    col_tot = np.where(fail, 0.0, A).sum(axis=0)
    surviving = {alleles[i] for i in range(1, n) if col_tot[i] > 0}
    return alleles[0], surviving


def stale_format_fields(fpath: str, keep=("GT", "AD")) -> list[str]:
    """FORMAT fields invalidated when a record's alleles are reshaped.

    Any per-allele / per-genotype field (Number A, R or G) other than GT and AD
    goes stale once harmonize edits the allele set — harmonize maintains AD but
    cannot recompute e.g. PL, so those must be dropped or a downstream
    ``bcftools merge`` fails with a FORMAT length mismatch. ``keep`` names the fields
    harmonize does maintain: every ``Number=R`` count (ADF, ADR) is padded like AD, so a
    caller that has not zeroed any depths can keep those too.
    """
    import pysam

    stale = []
    with pysam.VariantFile(fpath) as vcf:
        for fid, meta in vcf.header.formats.items():
            if fid in keep:
                continue
            if str(meta.number) in ("A", "R", "G"):
                stale.append(fid)
    return stale


def is_indel_context(rec) -> bool:
    """Whether a record is an indel-context record that does not belong in a SNP set.

    Catches records that ``bcftools view --exclude-type indels`` misses because
    their ALT is ``.`` (no alternate allele to type as an indel) yet they carry
    the ``INDEL`` INFO flag, plus any multi-base REF/ALT. Legitimate SNP records
    — including a monomorphic ``REF>.`` site with no INDEL flag — are kept, since
    those correctly get their ALT filled from the cross-file union.
    """
    if "INDEL" in rec.info:
        return True
    if len(rec.ref) > 1:
        return True
    for a in rec.alleles[1:]:
        if a != "." and len(a) != 1:
            return True
    return False


def _prefer(cand_n: int, prev_n: int) -> bool:
    """Duplicate-collapse rule: keep the record with the most real ALTs; ties
    keep the later record. This keeps the true SNP over an overlapping no-ALT /
    indel record emitted at the same start position."""
    return cand_n >= prev_n


def _pad_to_ref(alleles, ref: str, common_ref: str):
    """Re-express ``alleles`` (whose REF is ``ref``) against ``common_ref``.

    Two records at one position can have REFs of different length when one is an indel:
    ``A > T`` and ``ATT > A`` both sit at the same POS. Their alleles are only comparable
    once written against one REF, and the longer one is it -- ``A > T`` becomes
    ``ATT > TTT``, which is the same variant. ``bcftools merge`` does exactly this. Before
    this the union was ``[first file's REF] + sorted(every file's ALTs)`` and, with the SNP
    file first, came out as ``REF=A ALT=A,T``: an ALT equal to the reference, with the
    deletion carriers re-labelled as carrying it.

    ``common_ref`` must extend ``ref``; the caller has checked that. ``*`` and symbolic
    alleles are not sequences and are left alone.
    """
    pad = common_ref[len(ref):]
    if not pad:
        return list(alleles)
    return [a if (a == "*" or a.startswith("<") or a == ".") else a + pad for a in alleles]


def accumulate_union(files: list[str], min_ad: int, min_af: float, het_min_af: float,
                     drop_indels: bool = True, keep_ref_only: bool = False):
    """Pass 1: stream each file, clean, and collect real ALTs per site.

    Indel-context records are dropped by default (see :func:`is_indel_context`).
    ``keep_ref_only`` puts sites with no real ALT in any file into the union too, as
    ``[ref]``, so pass 2 writes them through as ``REF>.`` records instead of dropping
    them: a callset made over a list of positions is meant to answer at every position,
    and a reference call there is the answer.
    Within a file, records that share a ``(chrom, pos)`` are collapsed by keeping
    the one with the most real ALT alleles (see :func:`_prefer`) — so an
    overlapping no-ALT/indel record does not clobber the real SNP. Across files
    the surviving real ALTs are unioned.

    Returns ``(union, dup_positions, ambiguous, stats)``:
      * ``union``        — ``(chrom, pos) -> [ref, alt1, ...]`` for sites with a real ALT
      * ``dup_positions``— ``(file, chrom, pos)`` collapsed (SNP kept over no-ALT record)
      * ``ambiguous``    — ``(file, chrom, pos)`` where >1 record carried real ALTs
                           (genuinely un-normalized; needs `bcftools norm`)
      * ``stats``        — per-file cleaning counts plus the union tally, for the report:
        what cleaning actually did is the difference between a threshold that is doing
        useful work and one that is quietly discarding real alleles.
    """
    ref_of: dict = {}
    raw_of: dict = {}   # key -> [(ref, {alts}) per file]
    dup_positions: set = set()
    ambiguous: set = set()
    per_file_stats: dict = {}
    import pysam

    for fpath in files:
        per_file: dict = {}  # key -> (ref, {alts}, n_real)
        st = {"processed": 0, "indel_context": 0, "alts_removed": 0, "reduced_to_ref_only": 0}
        with pysam.VariantFile(fpath) as vcf:
            for rec in vcf:
                st["processed"] += 1
                if drop_indels and is_indel_context(rec):
                    st["indel_context"] += 1
                    continue
                # Pass 1 only needs the surviving-ALT set, so skip re-genotyping.
                ref, real = surviving_alleles(rec, min_ad, min_af)
                before = n_real_alts(rec)
                st["alts_removed"] += max(0, before - len(real))
                if before > 0 and not real:
                    st["reduced_to_ref_only"] += 1
                key = (rec.chrom, rec.pos)
                cand = (ref, real, len(real))
                if key in per_file:
                    dup_positions.add((fpath, rec.chrom, rec.pos))
                    if per_file[key][2] >= 1 and cand[2] >= 1:
                        ambiguous.add((fpath, rec.chrom, rec.pos))
                    if _prefer(cand[2], per_file[key][2]):
                        per_file[key] = cand
                else:
                    per_file[key] = cand
        st["sites"] = len(per_file)
        per_file_stats[fpath] = st
        for key, (ref, alts, _n) in per_file.items():
            # the longest REF at a position is the one every file's alleles get written
            # against; every other REF there has to be a prefix of it or the files do not
            # agree on the reference sequence, which no re-expression can fix
            have = ref_of.get(key)
            if have is None or len(ref) > len(have):
                ref_of[key] = ref
            raw_of.setdefault(key, []).append((ref, alts))

    ref_conflicts = []
    alts_of = {}
    for key, entries in raw_of.items():
        common = ref_of[key]
        merged = set()
        for ref, alts in entries:
            if not common.startswith(ref):
                ref_conflicts.append((key, ref, common))
                break
            merged.update(_pad_to_ref(alts, ref, common))
        else:
            alts_of[key] = merged
    if ref_conflicts:
        shown = "; ".join(f"{c}:{p} has REF {a!r} and {b!r}" for (c, p), a, b in ref_conflicts[:5])
        raise SystemExit(
            f"ERROR: {len(ref_conflicts)} position(s) have REF alleles that are not "
            f"prefixes of one another across the inputs ({shown}"
            f"{'; ...' if len(ref_conflicts) > 5 else ''}). The files do not agree on the "
            "reference sequence there, and no re-expression of the alleles can reconcile "
            "them. Check that every input was called against the same reference.")

    union = {k: [ref_of[k]] + sorted(alts) for k, alts in alts_of.items()
             if alts or keep_ref_only}
    with_alts = sum(1 for alts in alts_of.values() if alts)
    stats = {"per_file": per_file_stats, "union_sites": len(alts_of),
             "union_with_alts": with_alts, "union_dropped": len(alts_of) - len(union)}
    return union, dup_positions, ambiguous, stats


def harmonize_record_to_union(rec, union_alleles, het_min_af, out, *,
                              regenotype: bool = True, shapes: dict | None = None,
                              stale_info=STALE_INFO_FIELDS) -> bool:
    """Rewrite one cleaned record to the union allele set and write it.

    Every ``Number=R`` count field (AD, and ADF/ADR where present, in FORMAT and INFO) is
    re-laid-out on the union: a value the record already had moves to the allele's new
    slot, and an allele the record did not carry gets 0 -- the file was called over these
    reads and found none of it. ``regenotype`` decides what happens to GT: ``True`` (the
    cross-cohort default) re-calls from AD any sample at a record that gained alleles;
    ``False`` keeps the caller's genotype and only re-indexes it, which is right when the
    files are one caller's output split by sample rather than separate cohorts.

    Returns whether the record gained ALT alleles it did not carry, which is the
    interesting half of the tally: those are the records whose genotypes were recomputed.
    """
    if shapes is None:
        shapes = field_shapes(rec.header)
    stale = list(stale_info) + [f for f in shapes["info_ag"] if f not in stale_info]
    union_alts = union_alleles[1:]
    # the union may be written against a longer REF than this record's (see _pad_to_ref);
    # its ALTs only match the union's once written the same way
    current_alts = _pad_to_ref([a for a in rec.alleles[1:] if a != "."],
                               rec.ref, union_alleles[0])

    if current_alts == union_alts and rec.ref == union_alleles[0]:
        strip_stale_info(rec, stale)
        out.write(rec)
        return False

    current_alt_to_idx = {a: i + 1 for i, a in enumerate(current_alts)}
    union_to_current = [0 if i == 0 else current_alt_to_idx.get(a)
                        for i, a in enumerate(union_alleles)]
    old_to_new = {old: new for new, old in enumerate(union_to_current) if old is not None}

    # snapshot the per-allele arrays before touching rec.alleles (pysam reinterprets them)
    n_union = len(union_alleles)
    info_r, fmt_r = _snapshot_r(rec, shapes)
    gts = {sname: rec.samples[sname].get("GT", None) for sname in rec.samples}

    rec.alleles = tuple(union_alleles)
    alleles_added = any(idx is None for idx in union_to_current[1:])

    def relaid(old):
        return tuple(old[idx] if (idx is not None and idx < len(old)) else 0
                     for idx in union_to_current)

    for fid in shapes["info_r"]:
        if info_r[fid] is not None:
            rec.info[fid] = relaid(info_r[fid])
    for sname in rec.samples:
        sample = rec.samples[sname]
        for fid in shapes["format_r"]:
            if fid == "AD":
                continue
            v = fmt_r[sname][fid]
            sample[fid] = relaid(v) if v is not None else tuple([0] * n_union)
        old_ad = fmt_r[sname].get("AD")
        if old_ad is None:
            sample["AD"] = tuple([0] * n_union)
            sample["GT"] = (None, None)
            continue
        new_ad = relaid(old_ad)
        sample["AD"] = new_ad
        if alleles_added and regenotype:
            gt = regenotype_from_ad(list(new_ad), het_min_af)
            sample["GT"] = gt if gt is not None else (None, None)
        elif gts[sname] is not None:
            # same alleles in a different order, or the caller's call kept: re-index
            _remap_gt(sample, old_to_new)
    strip_stale_info(rec, stale)
    out.write(rec)
    return alleles_added


#: output-format code -> file extension
OUTPUT_EXT = {"v": ".vcf", "z": ".vcf.gz", "b": ".bcf"}


def harmonize_file(fpath: str, out_path: str, union: dict,
                   min_ad: int, min_af: float, het_min_af: float,
                   drop_indels: bool = True, *, regenotype: bool = True,
                   stale_info=STALE_INFO_FIELDS) -> dict:
    """Pass 2: stream a file, clean each record, and write it harmonized to union.

    Exactly one record is written per ``(chrom, pos)`` — the one with the most
    real ALT alleles (see :func:`_prefer`) — so duplicate input positions do not
    produce duplicate (merge-breaking) output positions and the real SNP is kept
    over an overlapping no-ALT/indel record. Inputs are assumed coordinate-sorted,
    so duplicates are consecutive.

    Always writes VCF text. Records whose allele count is reduced during cleaning
    must not be written straight to BCF: pysam does not shrink the ``Number=R``
    AD array to match, leaving a binary AD/allele mismatch that breaks downstream
    tools. Converting the VCF to BCF (e.g. via bcftools) regenerates AD cleanly.

    ``regenotype`` and ``stale_info`` are passed through to :func:`clean_record` and
    :func:`harmonize_record_to_union`.

    Returns ``{"written", "alts_added", "dropped_ref_only", "absent"}``. ``absent`` counts
    union sites this file holds no record for: harmonizing makes the files agree on *alleles*,
    not on which sites they contain, so those become whole-cohort missing genotypes after
    ``bcftools merge`` — see the note the command prints.
    """
    import pysam

    st = {"written": 0, "alts_added": 0, "dropped_ref_only": 0}
    seen: set = set()
    with pysam.VariantFile(fpath) as vcf:
        shapes = field_shapes(vcf.header)
        opts = dict(regenotype=regenotype, shapes=shapes, stale_info=stale_info)
        out = pysam.VariantFile(out_path, "w", header=vcf.header)
        held = None
        held_key = None
        held_n = -1
        for rec in vcf:
            if drop_indels and is_indel_context(rec):
                continue
            clean_record(rec, min_ad, min_af, het_min_af, **opts)
            key = (rec.chrom, rec.pos)
            if held is not None and key != held_key:
                _emit(held, held_key, union, het_min_af, out, st, seen, opts)
                held, held_n = None, -1
            if key in seen:
                # Pass 1 collapses duplicate positions over the whole file; this loop only
                # collapses *adjacent* ones. On a coordinate-sorted file the two agree. On
                # anything else they can pick different records, and the record emitted here
                # may name an allele the union was never told about -- which used to surface
                # as a bare KeyError deep in _remap_gt, or, with regenotype on, as a sample
                # silently re-called from an all-zero AD.
                raise SystemExit(
                    f"ERROR: {fpath} is not coordinate-sorted: {key[0]}:{key[1]} appears "
                    f"again after another position.\n"
                    "  harmonize builds its ALT union in one pass and rewrites in a second, "
                    "and the two only agree on sorted input.\n"
                    "  Run `bcftools sort` on it first."
                )
            cand_n = n_real_alts(rec)
            if held is None or _prefer(cand_n, held_n):
                held, held_key, held_n = rec, key, cand_n
        if held is not None:
            _emit(held, held_key, union, het_min_af, out, st, seen, opts)
        out.close()
    st["absent"] = len(union) - len(seen)
    return st


def _emit(rec, key, union, het_min_af, out, st: dict, seen: set, opts: dict) -> None:
    union_alleles = union.get(key)
    if union_alleles is None:
        st["dropped_ref_only"] += 1
        return  # site dropped (ref-only across all files)
    if harmonize_record_to_union(rec, union_alleles, het_min_af, out, **opts):
        st["alts_added"] += 1
    st["written"] += 1
    seen.add(key)
