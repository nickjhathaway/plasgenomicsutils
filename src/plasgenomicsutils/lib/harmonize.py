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


def strip_stale_info(rec) -> None:
    for field in STALE_INFO_FIELDS:
        try:
            del rec.info[field]
        except KeyError:
            pass


def clean_record(rec, min_ad: int, min_af: float, het_min_af: float) -> None:
    """Zero low-level ALT depths, drop empty ALTs, re-genotype (mutates ``rec``).

    Records reduced to ref-only are kept with ALT="." so the union step can
    still fill ALTs from other files. REF depth is never zeroed.
    """
    alleles = list(rec.alleles)
    n_alleles = len(alleles)
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
        new_alleles = [alleles[i] for i in range(n_alleles) if alts_to_keep[i]]
        if len(new_alleles) == 1:
            rec.alleles = (alleles[0], ".")
            for sname in rec.samples:
                ad = new_ad_per_sample.get(sname)
                ref_ad = ad[0] if ad is not None else 0
                rec.samples[sname]["AD"] = (ref_ad, 0)
                rec.samples[sname]["GT"] = (0, 0)
            strip_stale_info(rec)
            return
        n_new = len(new_alleles)
        rec.alleles = tuple(new_alleles)
        for sname in rec.samples:
            sample = rec.samples[sname]
            ad = new_ad_per_sample.get(sname)
            if ad is not None and len(ad) == n_alleles:
                new_ad = [ad[i] for i in range(n_alleles) if alts_to_keep[i]]
                sample["AD"] = tuple(new_ad)
                gt = regenotype_from_ad(new_ad, het_min_af)
                sample["GT"] = gt if gt is not None else (None, None)
            else:
                # Sample had no usable AD; still write a correct-length AD so the
                # record stays consistent (Number=R) after the allele count drops.
                sample["AD"] = tuple([0] * n_new)
                sample["GT"] = (None, None)
        strip_stale_info(rec)
    else:
        for sname in rec.samples:
            sample = rec.samples[sname]
            ad = new_ad_per_sample.get(sname)
            if ad is None or len(ad) != n_alleles:
                continue
            original_ad = list(sample.get("AD", []) or [])
            sample["AD"] = tuple(ad)
            if ad != original_ad:
                gt = regenotype_from_ad(ad, het_min_af)
                sample["GT"] = gt if gt is not None else (None, None)
        strip_stale_info(rec)


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


def stale_format_fields(fpath: str) -> list[str]:
    """FORMAT fields invalidated when a record's alleles are reshaped.

    Any per-allele / per-genotype field (Number A, R or G) other than GT and AD
    goes stale once harmonize edits the allele set — harmonize maintains AD but
    cannot recompute e.g. PL, so those must be dropped or a downstream
    ``bcftools merge`` fails with a FORMAT length mismatch.
    """
    import pysam

    stale = []
    with pysam.VariantFile(fpath) as vcf:
        for fid, meta in vcf.header.formats.items():
            if fid in ("GT", "AD"):
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


def accumulate_union(files: list[str], min_ad: int, min_af: float, het_min_af: float,
                     drop_indels: bool = True):
    """Pass 1: stream each file, clean, and collect real ALTs per site.

    Indel-context records are dropped by default (see :func:`is_indel_context`).
    Within a file, records that share a ``(chrom, pos)`` are collapsed by keeping
    the one with the most real ALT alleles (see :func:`_prefer`) — so an
    overlapping no-ALT/indel record does not clobber the real SNP. Across files
    the surviving real ALTs are unioned.

    Returns ``(union, dup_positions, ambiguous)``:
      * ``union``        — ``(chrom, pos) -> [ref, alt1, ...]`` for sites with a real ALT
      * ``dup_positions``— ``(file, chrom, pos)`` collapsed (SNP kept over no-ALT record)
      * ``ambiguous``    — ``(file, chrom, pos)`` where >1 record carried real ALTs
                           (genuinely un-normalized; needs `bcftools norm`)
    """
    ref_of: dict = {}
    alts_of: dict = {}
    dup_positions: set = set()
    ambiguous: set = set()
    import pysam

    for fpath in files:
        per_file: dict = {}  # key -> (ref, {alts}, n_real)
        with pysam.VariantFile(fpath) as vcf:
            for rec in vcf:
                if drop_indels and is_indel_context(rec):
                    continue
                # Pass 1 only needs the surviving-ALT set, so skip re-genotyping.
                ref, real = surviving_alleles(rec, min_ad, min_af)
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
        for key, (ref, alts, _n) in per_file.items():
            ref_of.setdefault(key, ref)
            alts_of.setdefault(key, set()).update(alts)

    union = {k: [ref_of[k]] + sorted(alts) for k, alts in alts_of.items() if alts}
    return union, dup_positions, ambiguous


def harmonize_record_to_union(rec, union_alleles, het_min_af, out) -> None:
    """Rewrite one cleaned record to the union allele set and write it."""
    union_alts = union_alleles[1:]
    current_alts = [a for a in rec.alleles[1:] if a != "."]

    if current_alts == union_alts:
        strip_stale_info(rec)
        out.write(rec)
        return

    current_alt_to_idx = {a: i + 1 for i, a in enumerate(current_alts)}
    union_to_current = [0 if i == 0 else current_alt_to_idx.get(a)
                        for i, a in enumerate(union_alleles)]

    # snapshot AD before touching rec.alleles (pysam reinterprets FORMAT on change)
    n_union = len(union_alleles)
    ad_snapshot = {}
    for sname in rec.samples:
        raw = rec.samples[sname].get("AD", None)
        ad_snapshot[sname] = None if (raw is None or any(v is None for v in raw)) else list(raw)

    rec.alleles = tuple(union_alleles)
    alleles_added = any(idx is None for idx in union_to_current[1:])

    for sname in rec.samples:
        sample = rec.samples[sname]
        old_ad = ad_snapshot[sname]
        if old_ad is None:
            sample["AD"] = tuple([0] * n_union)
            sample["GT"] = (None, None)
            continue
        new_ad = [old_ad[idx] if (idx is not None and idx < len(old_ad)) else 0
                  for idx in union_to_current]
        sample["AD"] = tuple(new_ad)
        if alleles_added:
            gt = regenotype_from_ad(new_ad, het_min_af)
            sample["GT"] = gt if gt is not None else (None, None)
    strip_stale_info(rec)
    out.write(rec)


#: output-format code -> file extension
OUTPUT_EXT = {"v": ".vcf", "z": ".vcf.gz", "b": ".bcf"}


def harmonize_file(fpath: str, out_path: str, union: dict,
                   min_ad: int, min_af: float, het_min_af: float,
                   drop_indels: bool = True) -> None:
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
    """
    import pysam

    with pysam.VariantFile(fpath) as vcf:
        out = pysam.VariantFile(out_path, "w", header=vcf.header)
        held = None
        held_key = None
        held_n = -1
        for rec in vcf:
            if drop_indels and is_indel_context(rec):
                continue
            clean_record(rec, min_ad, min_af, het_min_af)
            key = (rec.chrom, rec.pos)
            if held is not None and key != held_key:
                _emit(held, held_key, union, het_min_af, out)
                held, held_n = None, -1
            cand_n = n_real_alts(rec)
            if held is None or _prefer(cand_n, held_n):
                held, held_key, held_n = rec, key, cand_n
        if held is not None:
            _emit(held, held_key, union, het_min_af, out)
        out.close()


def _emit(rec, key, union, het_min_af, out) -> None:
    union_alleles = union.get(key)
    if union_alleles is None:
        return  # site dropped (ref-only across all files)
    harmonize_record_to_union(rec, union_alleles, het_min_af, out)
