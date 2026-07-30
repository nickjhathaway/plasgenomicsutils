"""Strand-bias sequencing-artifact diagnostics (SSE fake-het detection).

Illumina sequence-specific errors (SSE) at GC islands, G-homopolymers, and GGC
motifs cause lagging-strand dephasing that produces **fake heterozygous calls**
whose ALT reads sit almost entirely on one strand, at low base quality. In
CN-from-VAF analyses these land near the 0.25/0.33 fractions a real 3-4 copy locus
would show, so they mimic the very signal they corrupt and must be excluded. A real
het/CNV allele is strand-balanced; an artifact is strand-restricted.

This module provides, in rough order of pipeline usefulness:

  * :func:`strand_bias_verdict`   -- per-site/-sample call from ADF/ADR (VCF-only, portable)
  * :func:`scan_vcf_strand_bias`  -- batch the verdict over a biallelic VCF/BCF with ADF/ADR
  * :func:`per_read_records`      -- single-position read-level tabulation (needs the BAM)
  * :func:`summarize_reverse_alt` -- artifact summary of the reverse-strand ALT reads
  * :func:`format_read_check_report` -- human-readable report of the above
  * :func:`extract_alt_reads`     -- pull ALT-carrying reads for direct viewing (FASTQ + aligned FASTA)
  * :func:`scan_reference_triggers` -- flag SSE-prone motifs from the reference alone

Reference: Nakamura K, et al. "Sequence-specific error profile of Illumina
sequencers." Nucleic Acids Research 2011;39(13):e90.
"""

from __future__ import annotations

from collections import Counter
from math import log10
from statistics import median, quantiles

from .pileup import iter_site_reads, revcomp

# --------------------------------------------------------------------------- #
#  VCF-only path: per-strand VAF + Fisher strand bias from ADF/ADR             #
# --------------------------------------------------------------------------- #


def strand_bias_verdict(ref_fwd, alt_fwd, ref_rev, alt_rev, *,
                        min_minor_depth: int = 20, ratio_hard: float = 0.15,
                        ratio_soft: float = 0.30, sb_hard: float = 60.0,
                        alt_bq_median=None):
    """Judge one candidate het (one sample, one site) from strand-resolved depths.

    ``ADF`` gives ``(ref_fwd, alt_fwd)`` and ``ADR`` gives ``(ref_rev, alt_rev)``
    after ``bcftools norm -m-``. Returns a dict with per-strand VAF, the Fisher
    strand-bias Phred (``sb_phred``), the minor/major VAF ``ratio``, a boolean
    ``drop``, and the ``reasons`` it fired.

    The combined rule (conservative about false drops) drops the call if:

      * strand-restricted (``ratio < ratio_hard``) AND the minor strand had the
        depth to detect the alt (``min_minor_depth``) — so a genuinely shallow
        strand is not mistaken for bias, OR
      * extreme Fisher bias (``sb_phred > sb_hard``), OR
      * a softer strand skew (``ratio < ratio_soft``) corroborated by low alt-read
        base quality (``alt_bq_median < 15``), when that BQ is supplied.
    """
    from scipy.stats import fisher_exact

    dpf = ref_fwd + alt_fwd
    dpr = ref_rev + alt_rev
    vaf_f = alt_fwd / dpf if dpf else 0.0
    vaf_r = alt_rev / dpr if dpr else 0.0
    minor, major = sorted((vaf_f, vaf_r))
    ratio = (minor / major) if major > 0 else 1.0
    p = fisher_exact([[ref_fwd, alt_fwd], [ref_rev, alt_rev]])[1]
    sb_phred = -10 * log10(max(p, 1e-300))
    minor_depth = dpf if vaf_f <= vaf_r else dpr

    drop, reasons = False, []
    if ratio < ratio_hard and minor_depth >= min_minor_depth:
        drop = True
        reasons.append(f"strand-restricted (ratio={ratio:.2f}, minor_depth={minor_depth})")
    if sb_phred > sb_hard:
        drop = True
        reasons.append(f"fisher SB Phred={sb_phred:.0f}")
    if alt_bq_median is not None and ratio < ratio_soft and alt_bq_median < 15:
        drop = True
        reasons.append(f"low alt BQ={alt_bq_median} + strand skew (ratio={ratio:.2f})")
    return {"drop": drop, "vaf_fwd": vaf_f, "vaf_rev": vaf_r, "vaf": vaf_f, "ratio": ratio,
            "sb_phred": sb_phred, "minor_depth": minor_depth, "reasons": reasons}


def scan_vcf_strand_bias(input_vcf: str, *, min_vaf: float = 0.05,
                         min_minor_depth: int = 20, ratio_hard: float = 0.15,
                         ratio_soft: float = 0.30, sb_hard: float = 60.0):
    """Run :func:`strand_bias_verdict` over every candidate het in a biallelic VCF.

    The input must be **biallelic** (``bcftools norm -m-``) and carry per-strand
    ``FORMAT/ADF`` and ``FORMAT/ADR`` (from ``bcftools mpileup -a FORMAT/ADF,ADR``).
    A sample is a candidate at a site when its pooled alt VAF is >= ``min_vaf``.

    Yields one dict per candidate ``(site, sample)``: chrom, pos (1-based), ref, alt,
    sample, the four strand depths, the verdict fields, and ``reasons`` joined with
    ``;``. Read-level base quality is not available here, so the low-BQ rule never
    fires on this path (use the read-level tools to corroborate).
    """
    from cyvcf2 import VCF

    vcf = VCF(input_vcf)
    samples = vcf.samples
    try:
        for v in vcf:
            if len(v.ALT) != 1:
                raise ValueError(
                    f"{v.CHROM}:{v.POS} is not biallelic; run `bcftools norm -m-` first")
            try:
                adf = v.format("ADF")
                adr = v.format("ADR")
            except KeyError:  # cyvcf2 raises when the FORMAT field is not in the header
                adf = adr = None
            if adf is None or adr is None:
                raise ValueError(
                    "input lacks FORMAT/ADF and/or FORMAT/ADR; regenerate with "
                    "`bcftools mpileup -a FORMAT/ADF,FORMAT/ADR`")
            for i, sample in enumerate(samples):
                ref_fwd, alt_fwd = int(adf[i, 0]), int(adf[i, 1])
                ref_rev, alt_rev = int(adr[i, 0]), int(adr[i, 1])
                if min(ref_fwd, alt_fwd, ref_rev, alt_rev) < 0:  # missing strand depth
                    continue
                total = ref_fwd + alt_fwd + ref_rev + alt_rev
                alt_total = alt_fwd + alt_rev
                if total == 0 or alt_total == 0:
                    continue
                if alt_total / total < min_vaf:
                    continue
                verdict = strand_bias_verdict(
                    ref_fwd, alt_fwd, ref_rev, alt_rev, min_minor_depth=min_minor_depth,
                    ratio_hard=ratio_hard, ratio_soft=ratio_soft, sb_hard=sb_hard)
                yield {
                    "chrom": v.CHROM, "pos": v.POS, "ref": v.REF, "alt": v.ALT[0],
                    "sample": sample, "ref_fwd": ref_fwd, "alt_fwd": alt_fwd,
                    "ref_rev": ref_rev, "alt_rev": alt_rev,
                    "vaf_fwd": verdict["vaf_fwd"], "vaf_rev": verdict["vaf_rev"],
                    "vaf": alt_total / total, "ratio": verdict["ratio"],
                    "sb_phred": verdict["sb_phred"], "drop": verdict["drop"],
                    "reasons": ";".join(verdict["reasons"]),
                }
    finally:
        vcf.close()


# --------------------------------------------------------------------------- #
#  Read-level path: single position, single BAM                               #
# --------------------------------------------------------------------------- #

PER_READ_COLUMNS = ["read", "strand", "base", "is_alt", "bq", "mapq", "qpos",
                    "seq_cycle", "read_len", "dist_nearest_end", "sc_5p", "sc_3p",
                    "d_refstart", "d_refend", "has_sa", "seq_start", "seq_start_ext"]


def per_read_records(bam, chrom: str, pos0: int, ref_base: str, *, min_mapq: int = 0):
    """Tabulate, per read covering 0-based ``pos0``, the features that separate SSE
    artifacts from real alleles: strand, base, base quality, true sequencing cycle,
    distance to the nearest read end, soft-clip geometry, split-alignment (SA) tag,
    and the sequencing-start coordinate (a priming-boundary probe). Returns a list of
    dicts keyed by :data:`PER_READ_COLUMNS`.
    """
    ref_base = ref_base.upper()
    rows = []
    for pr, read, base, seq_cycle in iter_site_reads(bam, chrom, pos0, min_mapq=min_mapq):
        qpos = pr.query_position
        quals = read.query_qualities
        bq = quals[qpos] if quals is not None else -1
        length = len(read.query_sequence)
        rev = read.is_reverse
        dist_nearest_end = min(qpos, length - 1 - qpos)  # strand-agnostic read-position

        cig = read.cigartuples or []
        sc_5p = cig[0][1] if cig and cig[0][0] == 4 else 0
        sc_3p = cig[-1][1] if cig and cig[-1][0] == 4 else 0
        d_refstart = pos0 - read.reference_start
        d_refend = (read.reference_end - 1) - pos0
        has_sa = read.has_tag("SA")  # split alignment => chimera

        # Sequencing-start coordinate (where cycle 0 mapped). Reverse read: first
        # sequenced base is the HIGH-coord aligned end; forward read: the LOW-coord
        # end. The soft-clip-extended variant adds the clip on the sequencing-start
        # side, approximating the molecule/priming boundary.
        if rev:
            seq_start = read.reference_end
            seq_start_ext = seq_start + sc_3p
        else:
            seq_start = read.reference_start + 1
            seq_start_ext = seq_start - sc_5p

        rows.append(dict(
            read=read.query_name, strand="-" if rev else "+", base=base,
            is_alt=(base != ref_base), bq=bq, mapq=read.mapping_quality, qpos=qpos,
            seq_cycle=seq_cycle, read_len=length, dist_nearest_end=dist_nearest_end,
            sc_5p=sc_5p, sc_3p=sc_3p, d_refstart=d_refstart, d_refend=d_refend,
            has_sa=int(has_sa), seq_start=seq_start, seq_start_ext=seq_start_ext,
        ))
    return rows


def _pct(n, d):
    return 0.0 if d == 0 else 100.0 * n / d


def summarize_reverse_alt(rows, alt_base=None):
    """Summarize the reverse-strand ALT reads — where the artifact lives.

    Returns per-strand read/ALT counts and, for the reverse ALT reads, the base
    quality, sequencing-cycle spread, near-end fraction, soft-clip/SA fractions, and
    fragment-start (priming) spike stats. Pure over the :func:`per_read_records`
    output — no BAM access — so it is unit-testable on synthetic rows.
    """
    alt = alt_base.upper() if alt_base else None

    def is_alt(r):
        return r["is_alt"] and (alt is None or r["base"] == alt)

    fwd = [r for r in rows if r["strand"] == "+"]
    rev = [r for r in rows if r["strand"] == "-"]
    fwd_alt = [r for r in fwd if is_alt(r)]
    rev_alt = [r for r in rev if is_alt(r)]
    rev_ref = [r for r in rev if not r["is_alt"]]

    out = {
        "n_reads": len(rows), "n_fwd": len(fwd), "n_rev": len(rev),
        "n_fwd_alt": len(fwd_alt), "n_rev_alt": len(rev_alt),
        "vaf_fwd": _pct(len(fwd_alt), len(fwd)), "vaf_rev": _pct(len(rev_alt), len(rev)),
        "rev_alt": None,
    }
    if not rev_alt:
        return out

    bqs = [r["bq"] for r in rev_alt if r["bq"] >= 0]
    cyc = [r["seq_cycle"] for r in rev_alt]
    ends = [r["dist_nearest_end"] for r in rev_alt]
    read_len = median([r["read_len"] for r in rev_alt])
    cyc_iqr = None
    if len(cyc) > 3:
        q = quantiles(cyc, n=4)
        cyc_iqr = (q[0], q[2], q[2] - q[0])

    alt_starts = [r["seq_start"] for r in rev_alt]
    _, modal_alt, frac_alt = _coord_spike(alt_starts)
    _, modal_ext, frac_ext = _coord_spike([r["seq_start_ext"] for r in rev_alt])
    ref_modal = ref_frac = None
    if rev_ref:
        _, ref_modal, ref_frac = _coord_spike([r["seq_start"] for r in rev_ref])

    out["rev_alt"] = {
        "bq_median": (median(bqs) if bqs else None),
        "bq_min": (min(bqs) if bqs else None), "bq_max": (max(bqs) if bqs else None),
        "pct_bq_lt20": _pct(sum(b < 20 for b in bqs), len(bqs)) if bqs else 0.0,
        "cyc_median": median(cyc), "cyc_min": min(cyc), "cyc_max": max(cyc),
        "cyc_iqr": cyc_iqr, "read_len": read_len,
        "pct_near_end": _pct(sum(e < 5 for e in ends), len(ends)),
        "pct_softclip": _pct(sum(1 for r in rev_alt if r["sc_5p"] or r["sc_3p"]), len(rev_alt)),
        "pct_sa": _pct(sum(r["has_sa"] for r in rev_alt), len(rev_alt)),
        "modal_start": modal_alt, "modal_start_frac": frac_alt,
        "modal_start_ext": modal_ext, "modal_start_ext_frac": frac_ext,
        "ref_modal_start": ref_modal, "ref_modal_start_frac": ref_frac,
    }
    return out


def _coord_spike(coords):
    """Return ``(counter, modal_coord, modal_fraction)`` for a priming-spike check."""
    if not coords:
        return None, None, 0.0
    c = Counter(coords)
    modal, modal_n = c.most_common(1)[0]
    return c, modal, modal_n / len(coords)


def _texthist(vals, lo, hi, bins=20, width=40):
    if not vals:
        return "  (no values)"
    step = (hi - lo) / bins if hi > lo else 1
    counts = [0] * bins
    for v in vals:
        b = int((v - lo) / step) if step else 0
        counts[min(max(b, 0), bins - 1)] += 1
    mx = max(counts) or 1
    return "\n".join(f"  {lo + i * step:6.0f} |{'#' * int(round(width * c / mx))} {c}"
                     for i, c in enumerate(counts))


def format_read_check_report(rows, bam_name: str, chrom: str, pos1: int,
                             ref_base: str, alt_base=None) -> str:
    """Render the per-read tabulation as the diagnostic text report."""
    s = summarize_reverse_alt(rows, alt_base)
    alt = alt_base.upper() if alt_base else None
    lines = [
        f"=== {bam_name}  {chrom}:{pos1}  (ref={ref_base}"
        + (f", alt={alt}" if alt else "") + ") ===",
        f"  reads over site: {s['n_reads']}  (fwd {s['n_fwd']}, rev {s['n_rev']})",
        f"  fwd ALT reads          n={s['n_fwd_alt']}",
        f"  rev ALT reads          n={s['n_rev_alt']}",
        f"  strand VAF:  fwd {s['vaf_fwd']:.1f}%   rev {s['vaf_rev']:.1f}%",
    ]
    ra = s["rev_alt"]
    if ra is None:
        lines.append("  no reverse-strand ALT reads found at this position.")
        return "\n".join(lines)

    lines.append("\n  -- reverse-strand ALT reads (the artifact) --")
    if ra["bq_median"] is not None:
        lines.append(f"  base quality:  median {ra['bq_median']:.0f}  min {ra['bq_min']}  "
                     f"max {ra['bq_max']}  |  %BQ<20: {ra['pct_bq_lt20']:.0f}%")
    lines.append(f"  seq cycle:     median {ra['cyc_median']:.0f}  "
                 f"range [{ra['cyc_min']}, {ra['cyc_max']}]  (read len ~{ra['read_len']:.0f})")
    if ra["cyc_iqr"]:
        lo, hi, width = ra["cyc_iqr"]
        lines.append(f"                 IQR [{lo:.0f}, {hi:.0f}]  width {width:.0f} cycles")
    lines.append(f"  near read end (<5bp): {ra['pct_near_end']:.0f}%")
    lines.append(f"  soft-clipped reads:   {ra['pct_softclip']:.0f}%   "
                 f"SA-tag (chimera):  {ra['pct_sa']:.0f}%")
    lines.append("\n  seq-cycle distribution (reverse ALT reads):")
    cyc = [r["seq_cycle"] for r in rows
           if r["strand"] == "-" and r["is_alt"] and (alt is None or r["base"] == alt)]
    lines.append(_texthist(cyc, 0, max(ra["read_len"], max(cyc) + 1)))
    lines += [
        "\n  interpretation hints:",
        "   * cycle IQR narrow (<~15) & one cycle band  -> flow-cell / phasing / tile defect",
        "   * cycle spread across whole read, low BQ     -> DNA damage (strand-specific)",
        "   * high soft-clip% or SA%, POS near ref end   -> sWGA chimera junction",
        "\n  -- fragment-start (reverse reads: sequencing start = reference_end) --",
        f"  reverse ALT reads: modal start {ra['modal_start']} "
        f"({ra['modal_start_frac'] * 100:.0f}% of alt reads)   "
        f"soft-clip-extended modal {ra['modal_start_ext']} "
        f"({ra['modal_start_ext_frac'] * 100:.0f}%)",
    ]
    if ra["ref_modal_start"] is not None:
        lines.append(f"  reverse REF reads: modal start {ra['ref_modal_start']} "
                     f"({ra['ref_modal_start_frac'] * 100:.0f}% of ref reads)")
        if ra["modal_start_frac"] > 1.5 * ra["ref_modal_start_frac"]:
            lines.append("   * ALT reads MORE spike-concentrated than REF reads "
                         "-> artifact tied to reads from that priming site")
        else:
            lines.append("   * ALT and REF reads share start distribution "
                         "-> position-anchored (context) more than start-anchored")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
#  Read extraction for direct viewing                                         #
# --------------------------------------------------------------------------- #


def _qual_to_str(q):
    return "".join(chr(x + 33) for x in q) if q is not None else ""


def extract_alt_reads(bam, ref_fasta: str, chrom: str, pos1: int, out_prefix: str, *,
                      alt_base=None, window: int = 40, strand: str = "rev",
                      include_ref: int = 0, min_mapq: int = 0):
    """Write the reads carrying the ALT base at ``chrom:pos1`` for direct viewing.

    Produces two files, so an SSE artifact can be *seen*:

      * ``<prefix>.alt_reads.fastq`` — the ALT reads in original sequenced
        orientation (reverse reads reverse-complemented back), headers carrying
        strand / alt base / cycle / base quality / mate.
      * ``<prefix>.alt_reads.refaln.fasta`` — reference-anchored: every read padded
        over POS +/- ``window`` so the POS column lines up; the reference is the
        first record.

    ``strand`` selects ``rev`` (the usual artifact strand), ``fwd``, or ``both``.
    Up to ``include_ref`` reference-base reads (same strand) may be included for
    contrast. Returns ``(n_alt_reads, n_ref_reads, fastq_path, fasta_path)``; writes
    nothing and returns zero counts when no ALT reads are found.
    """
    import pysam

    pos0 = pos1 - 1
    fa = pysam.FastaFile(ref_fasta)
    try:
        ref_base = fa.fetch(chrom, pos0, pos0 + 1).upper()
        ref_win = fa.fetch(chrom, pos0 - window, pos0 + window + 1).upper()
    finally:
        fa.close()
    alt = alt_base.upper() if alt_base else None

    alt_reads, ref_reads = [], []
    for pr, read, base, seq_cycle in iter_site_reads(bam, chrom, pos0, min_mapq=min_mapq):
        rev = read.is_reverse
        if strand == "rev" and not rev:
            continue
        if strand == "fwd" and rev:
            continue
        quals = read.query_qualities
        bq = quals[pr.query_position] if quals is not None else -1
        rec = dict(read=read, base=base, cyc=seq_cycle, bq=bq, rev=rev)
        if base != ref_base and (alt is None or base == alt):
            alt_reads.append(rec)
        elif base == ref_base and len(ref_reads) < include_ref:
            ref_reads.append(rec)

    if not alt_reads:
        return 0, 0, None, None

    fastq_path = out_prefix + ".alt_reads.fastq"
    with open(fastq_path, "w") as fh:
        for rec in alt_reads:
            read = rec["read"]
            seq, q = read.query_sequence, read.query_qualities
            if rec["rev"]:
                seq = revcomp(seq)
                q = q[::-1] if q is not None else None
            mate = "/2" if read.is_read2 else ("/1" if read.is_read1 else "")
            hdr = (f"{read.query_name}{mate} strand={'-' if rec['rev'] else '+'} "
                   f"alt={rec['base']} cycle={rec['cyc']} bq={rec['bq']}")
            fh.write(f"@{hdr}\n{seq}\n+\n{_qual_to_str(q)}\n")

    win_start = pos0 - window
    win_len = 2 * window + 1

    def aligned_row(read):
        row = ["."] * win_len  # '.' = not covered
        for qpos, refpos in read.get_aligned_pairs():
            if refpos is None:  # insertion in read -> skip, keeping columns = reference
                continue
            col = refpos - win_start
            if 0 <= col < win_len:
                row[col] = "-" if qpos is None else read.query_sequence[qpos].upper()
        return "".join(row)

    fasta_path = out_prefix + ".alt_reads.refaln.fasta"
    with open(fasta_path, "w") as fh:
        fh.write(f">REFERENCE {chrom}:{pos1 - window}-{pos1 + window} "
                 f"(POS {pos1} is column {window}, ref {ref_base})\n{ref_win}\n")
        for tag, recs in (("alt", alt_reads), ("REF", ref_reads)):
            for rec in recs:
                read = rec["read"]
                mate = "/2" if read.is_read2 else ("/1" if read.is_read1 else "")
                label = f"{tag}{rec['base']}" if tag == "alt" else f"REF{rec['base']}"
                fh.write(f">{read.query_name}{mate}|{'-' if rec['rev'] else '+'}|"
                         f"{label}|cyc{rec['cyc']}|bq{rec['bq']}\n{aligned_row(read)}\n")

    return len(alt_reads), len(ref_reads), fastq_path, fasta_path


# --------------------------------------------------------------------------- #
#  Reference-context pre-flagging (cheap, deterministic)                      #
# --------------------------------------------------------------------------- #


def scan_reference_triggers(seq: str, start: int = 0, *, min_homopolymer: int = 3,
                            gc_window: int = 20, gc_threshold: float = 0.35):
    """Flag Nakamura SSE-trigger motifs in a reference window (no reads needed).

    SSE is sequence-driven, so risky positions can be flagged before looking at any
    read — and in the AT-rich *P. falciparum* genome (~19% GC) these motifs are rare
    and therefore specific. ``seq`` is the reference substring; ``start`` is its
    0-based genomic coordinate. Returns a list of ``{kind, start, end, detail}``
    features (coordinates genomic, 0-based half-open) for:

      * ``homopolymer`` — G or C runs >= ``min_homopolymer``,
      * ``motif`` — GGC / GCC,
      * ``gc_island`` — any ``gc_window``-bp window with GC fraction > ``gc_threshold``.
    """
    seq = seq.upper()
    n = len(seq)
    features = []

    i = 0
    while i < n:
        if seq[i] in "GC":
            j = i
            while j < n and seq[j] == seq[i]:
                j += 1
            if j - i >= min_homopolymer:
                features.append({"kind": "homopolymer", "start": start + i,
                                 "end": start + j, "detail": seq[i] * (j - i)})
            i = j
        else:
            i += 1

    for m in range(n - 2):
        tri = seq[m:m + 3]
        if tri in ("GGC", "GCC"):
            features.append({"kind": "motif", "start": start + m,
                             "end": start + m + 3, "detail": tri})

    if n >= gc_window:
        gc = [1 if b in "GC" else 0 for b in seq]
        run = sum(gc[:gc_window])
        for w in range(n - gc_window + 1):
            if w > 0:
                run += gc[w + gc_window - 1] - gc[w - 1]
            frac = run / gc_window
            if frac > gc_threshold:
                features.append({"kind": "gc_island", "start": start + w,
                                 "end": start + w + gc_window,
                                 "detail": f"GC={frac:.2f}"})
    return features
