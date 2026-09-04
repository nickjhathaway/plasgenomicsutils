"""How far apart the variants in a callset are, per chromosome.

Density says how many variants there are per unit of genome; spacing says how they are
arranged, and the two answer different questions. A chromosome with one variant every 300 bp
and one with dense clusters separated by deserts can report the same density, and only the
spread of the gaps tells them apart -- which is why this reports quartiles rather than a mean.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys

import numpy as np

from .bcftools import q, require

#: Percentiles reported for the gap distribution, in order.
PERCENTILES = (25, 50, 75, 95)

_REGION_RE = re.compile(r"^([^:]+)(?::(\d+)-(\d+))?$")


def parse_locus(spec: str) -> list[tuple[str, int, int]]:
    """``"chr:start-end,chr:start-end"`` -> 0-based half-open intervals.

    The region strings are **1-based inclusive**, as `bcftools -r` reads them and as anyone
    quoting a locus writes it. They are converted here, so everything downstream of this
    function is 0-based half-open like the rest of the package. A bare ``chr`` means the
    whole chromosome, which comes back as ``(chrom, 0, -1)`` -- the end filled in later from
    the callset's own header.
    """
    out = []
    for part in (p.strip() for p in spec.split(",")):
        if not part:
            continue
        m = _REGION_RE.match(part)
        if not m:
            # the comma separates regions, so it cannot also group digits -- worth saying,
            # because coordinates pasted out of a genome browser carry them
            hint = ("  Coordinates cannot contain commas: the comma separates regions.\n"
                    if re.search(r"\d$", part) or ":" not in part else "")
            raise SystemExit(f"ERROR: cannot read region {part!r}; expected chr or "
                             f"chr:start-end, 1-based as in bcftools -r.\n{hint}".rstrip())
        chrom, start, end = m.groups()
        if start is None:
            out.append((chrom, 0, -1))
            continue
        s, e = int(start), int(end)
        if e < s:
            raise SystemExit(f"ERROR: region {part!r} ends before it starts")
        out.append((chrom, s - 1, e))          # 1-based inclusive -> 0-based half-open
    if not out:
        raise SystemExit("ERROR: --locus was empty")
    return out


def read_bed(path: str) -> list[tuple[str, int, int]]:
    """A BED's intervals, already 0-based half-open."""
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 3:
                raise SystemExit(f"ERROR: {path} is not a BED (needs chrom, start, end)")
            out.append((f[0], int(f[1]), int(f[2])))
    if not out:
        raise SystemExit(f"ERROR: {path} has no intervals")
    return out


def merge_intervals(ivs: list[tuple[str, int, int]]) -> dict[str, list[tuple[int, int]]]:
    """Union the intervals per chromosome, so overlapping ones are not counted twice.

    ``--locus A:100-200,A:150-160`` asks about 101 bases, not 112: bcftools reports each
    variant once whatever the overlap, so the denominator has to agree with it or the
    density comes out low.
    """
    by: dict[str, list[tuple[int, int]]] = {}
    for chrom, s, e in ivs:
        by.setdefault(chrom, []).append((s, e))
    for chrom, spans in by.items():
        merged: list[tuple[int, int]] = []
        for s, e in sorted(spans):
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        by[chrom] = merged
    return by


def contig_lengths(path: str) -> dict[str, int]:
    """Contig lengths from the callset's own header."""
    require("bcftools")
    p = subprocess.run(f"bcftools view -h {q(path)}", shell=True, executable="/bin/bash",
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    out = {}
    for line in p.stdout.splitlines():
        if line.startswith("##contig="):
            m = re.search(r"ID=([^,>]+).*?length=(\d+)", line)
            if m:
                out[m.group(1)] = int(m.group(2))
    return out


def _is_indexed(path: str) -> bool:
    return any(os.path.exists(path + ext) for ext in (".csi", ".tbi"))


def _positions(path: str, locus: str | None, bed: str | None,
               quiet: bool = False) -> dict[str, np.ndarray]:
    """Variant positions per chromosome, 0-based, in file order.

    Regions go through ``-r``/``-R`` when the callset is indexed, which seeks straight to
    them, and through ``-t``/``-T`` when it is not. The targets form reads the whole file, so
    it is slower on a big callset -- but it is the difference between answering the question
    and refusing to, and plenty of working callsets are never indexed. Both forms read a
    region string as 1-based and a BED as 0-based, so the answer is the same either way.
    """
    require("bcftools")
    region = ""
    if locus or bed:
        indexed = _is_indexed(path)
        flag = ("-r" if indexed else "-t") if locus else ("-R" if indexed else "-T")
        region = f"{flag} {q(locus or bed)} "
        if not indexed and not quiet:
            sys.stderr.write(f"note: {path} has no index, so the regions are applied by "
                             f"scanning the whole file (bcftools {flag}). "
                             f"`bcftools index` makes this seek instead.\n")
    p = subprocess.run(f"bcftools query {region}-f '%CHROM\\t%POS0\\n' {q(path)}",
                       shell=True, executable="/bin/bash", stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise SystemExit(f"ERROR: bcftools query failed:\n{p.stderr.strip()}")
    by: dict[str, list[int]] = {}
    for line in p.stdout.splitlines():
        chrom, pos = line.split("\t")
        by.setdefault(chrom, []).append(int(pos))
    return {c: np.asarray(v, dtype=np.int64) for c, v in by.items()}


def variant_spacing(path: str, *, locus: str | None = None, bed: str | None = None,
                    bp_per_cm: float = 15000.0, quiet: bool = False) -> list[dict]:
    """Per-chromosome gap statistics and variant density, plus a pooled ``all`` row.

    Gaps are the distances between consecutive variants **within** a chromosome; the step
    across a chromosome boundary is not a gap and is never counted. A chromosome carrying one
    variant therefore has a count but no gap statistics, and the fields come back ``None``.

    ``locus`` is a bcftools region string (1-based, comma-separated); ``bed`` a BED file
    (0-based). The density denominator is the span actually asked about -- the union of those
    intervals per chromosome, or the contig's length from the header when neither is given.
    """
    if locus and bed:
        raise SystemExit("ERROR: pass --locus or --bed, not both")

    # regions are checked against the header before anything is read: on an unindexed callset
    # the query is a full scan, and a mistyped contig should not cost one to find out
    lengths = contig_lengths(path)
    if locus or bed:
        ivs = parse_locus(locus) if locus else read_bed(bed)
        unknown = sorted({c for c, _, _ in ivs} - set(lengths)) if lengths else []
        if unknown:
            raise SystemExit(
                f"ERROR: {', '.join(unknown)} is not a contig in {path}.\n"
                f"  It has: {', '.join(sorted(lengths)[:8])}"
                + (", ..." if len(lengths) > 8 else ""))
        ivs = [(c, s, lengths.get(c, 0) if e == -1 else e) for c, s, e in ivs]
        spans = {c: sum(e - s for s, e in v) for c, v in merge_intervals(ivs).items()}
    else:
        spans = dict(lengths)

    pos = _positions(path, locus, bed, quiet=quiet)

    # Rows follow what was *asked about*, not what was found: a region with no variants in it
    # is an answer -- often the interesting one -- and reporting nothing would read as a
    # failed query. Without regions there is nothing to ask about but the contigs with data.
    chroms = sorted(spans, key=_chrom_key) if (locus or bed) else sorted(pos, key=_chrom_key)
    rows, all_gaps, all_n, all_span = [], [], 0, 0
    for chrom in chroms:
        p = np.sort(pos.get(chrom, np.empty(0, dtype=np.int64)))
        gaps = np.diff(p)
        span = spans.get(chrom, 0)
        rows.append(_row(chrom, len(p), gaps, span, bp_per_cm))
        all_gaps.append(gaps)
        all_n += len(p)
        all_span += span
    if rows:
        rows.append(_row("all", all_n, np.concatenate(all_gaps) if all_gaps else np.array([]),
                         all_span, bp_per_cm))
    return rows


def _chrom_key(c: str):
    """Sort 1..14 numerically, then anything else (API, MIT) alphabetically after."""
    m = re.search(r"(\d+)", c)
    return (0, int(m.group(1)), "") if m else (1, 0, c)


def _row(chrom: str, n: int, gaps: np.ndarray, span: int, bp_per_cm: float) -> dict:
    row = {"chrom": chrom, "n_variants": n, "n_gaps": int(gaps.size), "span_bp": int(span),
           "variants_per_cm": round(n / (span / bp_per_cm), 3) if span else None,
           "min": None, "max": None}
    for pc in PERCENTILES:
        row[f"p{pc}"] = None
    if gaps.size:
        row["min"] = int(gaps.min())
        row["max"] = int(gaps.max())
        for pc, v in zip(PERCENTILES, np.percentile(gaps, PERCENTILES)):
            row[f"p{pc}"] = int(round(v))
    return row


#: Column order of the emitted table.
COLUMNS = ("chrom", "n_variants", "n_gaps", "span_bp", "variants_per_cm",
           "min", "p25", "p50", "p75", "p95", "max")


def write_table(rows: list[dict], out: str | None = None) -> None:
    """Write the rows as TSV, to ``out`` or stdout."""
    import sys

    fh = open(out, "w") if out else sys.stdout
    try:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in COLUMNS) + "\n")
    finally:
        if out:
            fh.close()
