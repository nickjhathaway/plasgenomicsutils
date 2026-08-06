"""Per-sample sequencing depth from BAM/CRAM, and cross-sample coverage dropouts.

Two questions this answers, both for QC:

  * **Is this sample deep enough?** Per chromosome and genome-wide: mean, median, SD and
    quartiles of depth, plus the *breadth* at a set of thresholds (the fraction of bases
    at >= 1x, 5x, 10x, 20x ...). Optionally restricted to a BED, which is the usual way
    to ask the question -- subtelomeric and hypervariable regions drag every statistic
    down, so the honest denominator is the core genome.
  * **Is this *region* ever covered, in anyone?** Selective whole-genome amplification
    (sWGA) does not amplify uniformly, and a region that no sample amplifies looks
    exactly like a region with no variation. :func:`dropout_regions` finds windows that
    are below depth in nearly every sample, so they can be excluded rather than silently
    read as invariant.

Everything is derived from one intermediate: a **depth histogram** (``hist[d]`` = number
of bases at depth ``d``). A histogram is enough for the mean, SD and any quantile or
threshold, and it stays small no matter how deep the sample, so a whole cohort can be
summarised without ever holding per-base depth for more than one chunk at a time.

Depth itself comes from one of two engines. **They do not define depth the same way**,
so the engine is recorded in the output and should be held fixed across a cohort:

  * ``pysam`` -- :meth:`pysam.AlignmentFile.count_coverage`. Counts *reads*, matching
    ``samtools depth`` base for base (verified on real WGS). Where the two mates of a
    fragment overlap, both are counted.
  * ``mosdepth`` -- much faster on whole-genome BAMs, and counts *fragments*: an
    overlapping mate pair contributes 1, not 2. On real Pf WGS this runs 2-3% below the
    read-level figure, everywhere, never above it.

Neither is wrong. Fragment depth is the better measure of independent evidence, since
two mates of one molecule are one observation; read depth is what most tools report.
What matters is not mixing them, so ``auto`` resolves once for the whole run and the
chosen engine is written into every row.

Coordinates are 0-based half-open throughout, matching the rest of the package (see
:mod:`plasgenomicsutils.lib.intervals`) and BED itself.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Read flags skipped when counting depth: unmapped, secondary, QC-fail, duplicate.
# Matches `samtools depth` / pysam's read_callback="all" default.
_SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400

DEFAULT_THRESHOLDS = (1, 5, 10, 20)
DEFAULT_WINDOW = 1000
DEFAULT_CHUNK = 5_000_000

#: A window counts as "uncovered" in a sample below this mean depth.
DROPOUT_MIN_DEPTH = 5.0
#: ...and the window is flagged when at least this fraction of samples are uncovered.
DROPOUT_MIN_FRAC_SAMPLES = 0.9


# --------------------------------------------------------------------------- #
#  Region sets                                                                 #
# --------------------------------------------------------------------------- #


def load_bed(path):
    """Read a BED into ``{chrom: [(start, end), ...]}``, merged and sorted.

    Overlapping and adjacent intervals are merged, so every base is counted once no
    matter how the file was written. Returns ``{}`` when ``path`` is falsy.
    """
    if not path:
        return {}
    raw: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 3:
                raise ValueError(f"BED line needs at least 3 columns: {line.rstrip()}")
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return {c: _merge(iv) for c, iv in raw.items()}


def _merge(intervals):
    out: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out


def _region_mask(intervals, chunk_start, chunk_end):
    """Boolean mask over ``[chunk_start, chunk_end)``; ``None`` intervals means all bases."""
    n = chunk_end - chunk_start
    if intervals is None:
        return None
    mask = np.zeros(n, dtype=bool)
    for start, end in intervals:
        if end <= chunk_start or start >= chunk_end:
            continue
        mask[max(start, chunk_start) - chunk_start:min(end, chunk_end) - chunk_start] = True
    return mask


# --------------------------------------------------------------------------- #
#  Depth accumulation                                                          #
# --------------------------------------------------------------------------- #


@dataclass
class ChromDepth:
    """Depth for one chromosome: a histogram plus per-window sums."""

    chrom: str
    length: int
    window: int
    hist: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.int64))
    win_depth_sum: np.ndarray = None      # sum of depth over the accessible bases
    win_bases: np.ndarray = None          # accessible bases per window

    def __post_init__(self):
        n_win = max(1, -(-self.length // self.window))
        if self.win_depth_sum is None:
            self.win_depth_sum = np.zeros(n_win, dtype=np.float64)
        if self.win_bases is None:
            self.win_bases = np.zeros(n_win, dtype=np.int64)

    def add(self, depth, offset, mask=None):
        """Fold a per-base depth array starting at ``offset`` into the accumulators."""
        if mask is not None:
            idx = np.flatnonzero(mask)
            if not idx.size:
                return
            depth = depth[idx]
            pos = idx + offset
        else:
            pos = np.arange(offset, offset + depth.size)
        counts = np.bincount(depth)
        if counts.size > self.hist.size:
            self.hist = np.pad(self.hist, (0, counts.size - self.hist.size))
        self.hist[:counts.size] += counts
        w = pos // self.window
        n_win = self.win_bases.size
        self.win_depth_sum += np.bincount(w, weights=depth, minlength=n_win)[:n_win]
        self.win_bases += np.bincount(w, minlength=n_win)[:n_win]


def summarise_histogram(hist, thresholds=DEFAULT_THRESHOLDS):
    """Mean, SD, quartiles and breadth-at-threshold from a depth histogram.

    ``hist[d]`` is the number of bases at depth ``d``. Quantiles are the smallest depth
    whose cumulative base count reaches the quantile, i.e. the usual discrete definition
    (no interpolation), so the median of an all-zero region is 0 rather than an average.
    """
    hist = np.asarray(hist, dtype=np.int64)
    n = int(hist.sum())
    out = {"bases": n, "mean": np.nan, "sd": np.nan, "min": np.nan, "median": np.nan,
           "q1": np.nan, "q3": np.nan, "max": np.nan, "pct_zero": np.nan}
    for t in thresholds:
        out[f"pct_ge_{t}x"] = np.nan
    if n == 0:
        return out

    depths = np.arange(hist.size, dtype=np.float64)
    mean = float((depths * hist).sum() / n)
    var = float((depths**2 * hist).sum() / n) - mean**2
    cum = np.cumsum(hist)

    def q(p):
        return float(np.searchsorted(cum, p * n, side="left"))

    nonzero = np.flatnonzero(hist)
    out.update(mean=mean, sd=float(np.sqrt(max(var, 0.0))),
               min=float(nonzero[0]), max=float(nonzero[-1]),
               median=q(0.5), q1=q(0.25), q3=q(0.75),
               pct_zero=100.0 * float(hist[0]) / n)
    # breadth: bases at >= t, straight off the tail of the cumulative count
    for t in thresholds:
        at_least = n - (int(cum[t - 1]) if t - 1 < cum.size else n)
        out[f"pct_ge_{t}x"] = 100.0 * at_least / n
    return out


# --------------------------------------------------------------------------- #
#  Engines                                                                     #
# --------------------------------------------------------------------------- #


def mosdepth_available():
    """Path to a usable ``mosdepth``, or ``None``."""
    return shutil.which("mosdepth")


def resolve_engine(engine, min_baseq=0):
    """Turn ``auto`` into the engine that will actually run.

    mosdepth has no base-quality filter, so a ``min_baseq`` request keeps the work on
    pysam (and is an error if mosdepth was named explicitly).
    """
    if engine not in ("auto", "pysam", "mosdepth"):
        raise ValueError(f"unknown engine: {engine}")
    if engine == "mosdepth":
        if not mosdepth_available():
            raise SystemExit("--engine mosdepth: the mosdepth binary is not on PATH")
        if min_baseq:
            raise SystemExit("--engine mosdepth cannot apply --min-baseq; use --engine pysam")
        return "mosdepth"
    if engine == "auto":
        return "mosdepth" if (mosdepth_available() and not min_baseq) else "pysam"
    return engine


def _chrom_lengths(bam_path, reference=None):
    import pysam

    with pysam.AlignmentFile(bam_path, reference_filename=reference) as bam:
        return dict(zip(bam.references, bam.lengths))


def _depth_pysam(bam_path, chrom, start, end, min_mapq, min_baseq, reference=None):
    """Per-base depth over ``[start, end)`` via pysam."""
    import pysam

    with pysam.AlignmentFile(bam_path, reference_filename=reference) as bam:
        if min_mapq > 0:
            # a Python callback runs per read, so this is the slow path -- only taken
            # when a MAPQ floor is actually asked for
            def keep(read):
                return not (read.flag & _SKIP_FLAGS) and read.mapping_quality >= min_mapq

            counts = bam.count_coverage(chrom, start, end, quality_threshold=min_baseq,
                                        read_callback=keep)
        else:
            counts = bam.count_coverage(chrom, start, end, quality_threshold=min_baseq,
                                        read_callback="all")
    return np.sum(np.asarray(counts, dtype=np.int64), axis=0)


def _run_mosdepth(bam_path, out_prefix, threads=1, min_mapq=0, reference=None):
    cmd = [mosdepth_available(), "-t", str(threads)]
    if min_mapq:
        cmd += ["--mapq", str(min_mapq)]
    if reference:
        cmd += ["--fasta", reference]
    cmd += [str(out_prefix), str(bam_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"mosdepth failed ({proc.returncode}):\n{proc.stderr.strip()}")
    return Path(f"{out_prefix}.per-base.bed.gz")


def _depth_from_per_base(path, lengths):
    """Expand mosdepth's run-length-encoded per-base BED into per-chromosome arrays."""
    per_chrom: dict[str, np.ndarray] = {}
    runs: dict[str, list[tuple[int, int, int]]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, d = line.split("\t")
            runs.setdefault(c, []).append((int(s), int(e), int(d)))
    for c, rl in runs.items():
        s = np.fromiter((r[0] for r in rl), dtype=np.int64, count=len(rl))
        e = np.fromiter((r[1] for r in rl), dtype=np.int64, count=len(rl))
        d = np.fromiter((r[2] for r in rl), dtype=np.int64, count=len(rl))
        arr = np.zeros(lengths.get(c, int(e[-1])), dtype=np.int64)
        if s[0] == 0 and np.array_equal(s[1:], e[:-1]):
            arr[:int(e[-1])] = np.repeat(d, e - s)      # contiguous: rebuild in one go
        else:
            for start, end, depth in zip(s, e, d):      # gapped: fill run by run
                arr[start:end] = depth
        per_chrom[c] = arr
    return per_chrom


# --------------------------------------------------------------------------- #
#  Per-sample coverage                                                         #
# --------------------------------------------------------------------------- #


def sample_coverage(bam_path, sample=None, regions=None, window=DEFAULT_WINDOW,
                    thresholds=DEFAULT_THRESHOLDS, chroms=None, engine="auto",
                    min_mapq=0, min_baseq=0, threads=1, reference=None,
                    chunk=DEFAULT_CHUNK):
    """Depth statistics for one BAM/CRAM.

    Args:
        bam_path: Indexed BAM/CRAM.
        sample: Name for the output rows (default: the file stem).
        regions: ``{chrom: [(start, end), ...]}`` from :func:`load_bed`, or ``None``
            for the whole genome. Statistics count only bases inside the region set,
            so a core-genome BED gives core-genome depth.
        window: Window size in bp for the per-window means used by
            :func:`dropout_regions`.
        thresholds: Depth cutoffs to report breadth for.
        chroms: Restrict to these contigs (default: all in the header, or all named
            in ``regions`` when one is given).
        engine: ``auto`` (mosdepth when installed), ``pysam`` or ``mosdepth``. These
            count different things -- see the module docstring -- so the resolved engine
            is returned in every row and should be held fixed across a cohort.

    Returns:
        ``(per_chrom, windows)``: a list of per-chromosome stat dicts with a final
        ``genome`` row, and a list of per-window dicts. Windows with no accessible
        bases are dropped.
    """
    bam_path = str(bam_path)
    sample = sample or Path(bam_path).name.split(".")[0]
    engine = resolve_engine(engine)
    lengths = _chrom_lengths(bam_path, reference)
    if chroms is None:
        chroms = [c for c in lengths if (not regions or c in regions)]
    else:
        unknown = [c for c in chroms if c not in lengths]
        if unknown:
            raise SystemExit(f"{bam_path}: contig(s) not in the header: {', '.join(unknown)}")

    acc: dict[str, ChromDepth] = {}
    if engine == "mosdepth":
        with tempfile.TemporaryDirectory() as td:
            per_base = _run_mosdepth(bam_path, Path(td) / sample, threads=threads,
                                     min_mapq=min_mapq, reference=reference)
            depth_by_chrom = _depth_from_per_base(per_base, lengths)
        for c in chroms:
            depth = depth_by_chrom.get(c)
            if depth is None:
                depth = np.zeros(lengths[c], dtype=np.int64)
            cd = ChromDepth(c, lengths[c], window)
            cd.add(depth, 0, _region_mask(regions.get(c, []) if regions else None,
                                          0, depth.size))
            acc[c] = cd
    else:
        for c in chroms:
            cd = ChromDepth(c, lengths[c], window)
            for start in range(0, lengths[c], chunk):
                end = min(start + chunk, lengths[c])
                depth = _depth_pysam(bam_path, c, start, end, min_mapq, min_baseq,
                                     reference)
                cd.add(depth, start,
                       _region_mask(regions.get(c, []) if regions else None, start, end))
            acc[c] = cd

    per_chrom = []
    genome_hist = np.zeros(1, dtype=np.int64)
    for c in chroms:
        cd = acc[c]
        row = {"sample": sample, "chrom": c, "engine": engine}
        row.update(summarise_histogram(cd.hist, thresholds))
        per_chrom.append(row)
        if cd.hist.size > genome_hist.size:
            genome_hist = np.pad(genome_hist, (0, cd.hist.size - genome_hist.size))
        genome_hist[:cd.hist.size] += cd.hist
    row = {"sample": sample, "chrom": "genome", "engine": engine}
    row.update(summarise_histogram(genome_hist, thresholds))
    per_chrom.append(row)

    windows = []
    for c in chroms:
        cd = acc[c]
        keep = np.flatnonzero(cd.win_bases > 0)
        for i in keep:
            windows.append({
                "sample": sample, "chrom": c,
                "start": int(i) * window,
                "end": min((int(i) + 1) * window, cd.length),
                "bases": int(cd.win_bases[i]),
                "mean_depth": float(cd.win_depth_sum[i] / cd.win_bases[i]),
            })
    return per_chrom, windows


# --------------------------------------------------------------------------- #
#  Cross-sample dropouts                                                       #
# --------------------------------------------------------------------------- #


def dropout_regions(windows, min_depth=DROPOUT_MIN_DEPTH,
                    min_frac_samples=DROPOUT_MIN_FRAC_SAMPLES, merge_gap=0,
                    min_length=0):
    """Windows that almost no sample covers, merged into regions.

    Args:
        windows: DataFrame with ``sample, chrom, start, end, mean_depth`` -- the second
            return value of :func:`sample_coverage`, concatenated over samples.
        min_depth: A sample is "uncovered" in a window below this mean depth.
        min_frac_samples: Flag a window when at least this fraction of samples are
            uncovered in it.
        merge_gap: Join flagged regions separated by at most this many bp.
        min_length: Drop merged regions shorter than this.

    Returns:
        A DataFrame of merged regions with the fraction of samples uncovered and the
        depth across samples, worst (most widely uncovered) first.
    """
    import pandas as pd

    if windows.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "length", "n_windows",
                                     "frac_samples_uncovered", "n_samples",
                                     "median_depth_across_samples",
                                     "mean_depth_across_samples"])
    w = windows.copy()
    w["uncovered"] = w["mean_depth"] < min_depth
    per_window = w.groupby(["chrom", "start", "end"], as_index=False).agg(
        n_samples=("sample", "nunique"),
        frac_samples_uncovered=("uncovered", "mean"),
        median_depth_across_samples=("mean_depth", "median"),
        mean_depth_across_samples=("mean_depth", "mean"),
    )
    flagged = per_window[per_window["frac_samples_uncovered"] >= min_frac_samples]
    if flagged.empty:
        return flagged.assign(length=0, n_windows=0).iloc[0:0]

    rows = []
    for chrom, sub in flagged.sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        block = None
        for r in sub.itertuples(index=False):
            if block is not None and r.start - block["end"] <= merge_gap:
                block["end"] = max(block["end"], r.end)
                block["parts"].append(r)
            else:
                if block is not None:
                    rows.append(block)
                block = {"chrom": chrom, "start": r.start, "end": r.end, "parts": [r]}
        if block is not None:
            rows.append(block)

    out = []
    for b in rows:
        parts = b["parts"]
        length = b["end"] - b["start"]
        if length < min_length:
            continue
        out.append({
            "chrom": b["chrom"], "start": b["start"], "end": b["end"], "length": length,
            "n_windows": len(parts),
            "n_samples": int(max(p.n_samples for p in parts)),
            "frac_samples_uncovered": float(
                np.mean([p.frac_samples_uncovered for p in parts])),
            "median_depth_across_samples": float(
                np.median([p.median_depth_across_samples for p in parts])),
            "mean_depth_across_samples": float(
                np.mean([p.mean_depth_across_samples for p in parts])),
        })
    df = pd.DataFrame(out)
    if df.empty:
        return df
    return df.sort_values(["frac_samples_uncovered", "length"],
                          ascending=[False, False]).reset_index(drop=True)


def annotate_regions(regions, genes):
    """Add a comma-separated ``genes`` column naming the genes each region overlaps.

    Args:
        regions: DataFrame with ``chrom, start, end``.
        genes: DataFrame with ``name, chr, start, end`` (0-based half-open).
    """
    import pandas as pd

    if regions.empty:
        return regions.assign(genes="")
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for g in genes.itertuples(index=False):
        by_chrom.setdefault(getattr(g, "chr"), []).append(
            (int(g.start), int(g.end), str(g.name)))
    names = []
    for r in regions.itertuples(index=False):
        hits = [n for s, e, n in by_chrom.get(r.chrom, []) if s < r.end and e > r.start]
        names.append(",".join(sorted(set(hits))))
    return regions.assign(genes=pd.Series(names, index=regions.index))
