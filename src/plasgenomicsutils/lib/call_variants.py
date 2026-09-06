"""Call variants with `bcftools mpileup | bcftools call`, annotated for the QC filter.

Two things this exists for.

**The annotations have to match what the filter reads.** `hard_qc_filter --caller bcftools`
tests `FS`, `RPBZ`, `SCBZ`, `MQBZ` and `MQSBZ`; only some of those come out of `mpileup` by
default, and a comparison against a tag that is not there is simply false -- so a callset
made without them passes the filter untouched and says nothing. Calling through here uses
the same annotation list the filter checks for, so the two cannot drift apart.

**Calling a region list is embarrassingly parallel, and `--threads` does not do it.**
`bcftools mpileup --threads` only parallelises compression; the pileup itself is one core.
The way to use a machine is to split the region list and run a job per chunk, which is what
this does: split, call each chunk concurrently, then `bcftools concat` the parts back into
one file. With no region list there is nothing to split and it runs a single job.

Splitting is not quite bit-for-bit, and it is worth knowing exactly how. The same positions
come out, with the same genotypes, depths, allele depths and bias statistics; **QUAL can
shift by a few points at a handful of records** -- indels and the odd SNP beside one --
because mpileup computes indel likelihoods and BAQ from the reads around a position, and
which neighbours share a chunk changes with the split. This is `bcftools mpileup -R`
itself, not this wrapper: cutting a region file in half by hand and concatenating the two
calls reproduces it. Nothing downstream in this package reads QUAL, but a QUAL cutoff of
your own is the one thing that could move.

Splitting a region file has one trap worth knowing about: **bcftools reads the coordinate
convention off the file extension.** A ``.bed`` is 0-based half-open, anything else is
1-based ``CHROM POS``. Chunks are therefore written with the same extension as the file
they came from, or the same list would mean two different things depending on how it was
split.

One BAM per sample is the usual arrangement here, so ``--ignore-RG`` is the default. That
names each sample after the *path* it was called from, which nothing downstream wants, so
the samples are renamed to the file name with ``sample_suffix`` removed as the last step.
"""

from __future__ import annotations

import math
import os
import pathlib
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

from .bcftools import index_vcf, out_flag, q, require, sample_names, sh
from .vcf_filters import BCFTOOLS_MPILEUP_ANNOTATIONS


def _region_lines(path: str) -> list[str]:
    """Region records, ignoring blanks and comments (which bcftools also skips)."""
    with open(path) as fh:
        return [ln for ln in fh
                if ln.strip() and not ln.lstrip().startswith(("#", "track", "browser"))]


def split_regions(path: str, outdir: str, *, n_chunks: int | None = None,
                  chunk_size: int | None = None) -> list[str]:
    """Split a region file into pieces, keeping its extension.

    ``chunk_size`` wins where given; otherwise the file is divided into ``n_chunks`` as
    evenly as the line count allows. Returns the chunk paths, and an empty list when the
    file holds a single region (nothing to gain by splitting one job into one job).
    """
    lines = _region_lines(path)
    if not lines:
        raise SystemExit(f"call_variants: {path} lists no regions")
    if chunk_size is None:
        n = max(1, int(n_chunks or 1))
        chunk_size = math.ceil(len(lines) / n)
    chunk_size = max(1, int(chunk_size))
    if chunk_size >= len(lines):
        return []

    # bcftools decides 0-based BED vs 1-based CHROM/POS from the extension, so a chunk of
    # a .bed has to stay a .bed
    base = os.path.basename(path)
    stem, ext = os.path.splitext(base)
    if ext == ".gz":                                   # .bed.gz -> keep both suffixes
        stem, inner = os.path.splitext(stem)
        ext = inner + ext
    out = []
    for i in range(0, len(lines), chunk_size):
        p = os.path.join(outdir, f"{stem}.chunk{i // chunk_size:04d}{ext}")
        with open(p, "w") as fh:
            fh.writelines(lines[i:i + chunk_size])
        out.append(p)
    return out


def sample_name_for(path: str, suffix: str | None = ".bam") -> str:
    """The sample name a BAM should carry: its file name, without ``suffix``.

    ``bcftools mpileup --ignore-RG`` names each sample after the path it was given, so a
    callset comes out with names like ``/tank/.../s1.sorted.dup.pf.bam``. Stripping the
    directory and a known suffix gives back what the sample is actually called; pass the
    whole trailing part you want gone (``--sample-suffix .sorted.dup.pf.bam``), since only
    you know which of it is pipeline bookkeeping and which is the name.
    """
    name = os.path.basename(str(path).strip())
    if suffix and name.endswith(suffix):
        name = name[: -len(suffix)]
    return name


def sample_rename_map(names: list[str], suffix: str | None = ".bam") -> dict[str, str]:
    """Old name -> new name, refusing to produce a collision or an empty name."""
    out = {}
    for old in names:
        new = sample_name_for(old, suffix)
        if not new:
            raise SystemExit(
                f"call_variants: stripping {suffix!r} from {old!r} leaves nothing to call "
                "the sample; pass a different --sample-suffix")
        out[old] = new
    clashes = {n for n in out.values() if list(out.values()).count(n) > 1}
    if clashes:
        offenders = sorted(o for o, n in out.items() if n in clashes)
        raise SystemExit(
            "call_variants: these alignments would end up with the same sample name after "
            f"stripping {suffix!r}: " + ", ".join(offenders) +
            "\n  Pass a --sample-suffix that keeps them apart, or --no-rename-samples to "
            "leave the names as bcftools wrote them.")
    return out


def _rename_samples(path: str, suffix: str | None) -> dict[str, str]:
    """Rename in place via `bcftools reheader`, returning the mapping applied."""
    mapping = sample_rename_map(sample_names(path), suffix)
    if all(o == n for o, n in mapping.items()):
        return {}
    fd, tmp_map = tempfile.mkstemp(suffix=".samples.txt")
    with os.fdopen(fd, "w") as fh:
        for old, new in mapping.items():
            fh.write(f"{old}\t{new}\n")
    tmp_out = f"{path}.reheader{os.path.splitext(path)[1]}"
    try:
        sh(f"bcftools reheader --samples {q(tmp_map)} {q(path)} -o {q(tmp_out)}",
           tools=("bcftools",))
        os.replace(tmp_out, path)
        index_vcf(path)
    finally:
        os.unlink(tmp_map)
        if os.path.exists(tmp_out):
            os.unlink(tmp_out)
    return mapping


def _rename_cmds(ignore_rg: bool, suffix: str | None, out: str) -> list[str]:
    """The reheader step, for --dry-run; the mapping needs the called file to exist."""
    if not (ignore_rg and suffix is not None):
        return []
    return [f"bcftools reheader --samples <path -> name, stripping {suffix!r}> "
            f"{q(out)} -o {q(out)}  # applied in place"]


def _mpileup_call_cmd(*, ref: str, bams: list[str] | None, bam_list: str | None,
                      out: str, regions: str | None, annotations: str, ploidy: str,
                      ignore_rg: bool, skip_indels: bool, max_depth: int | None,
                      min_mapq: int | None, min_baseq: int | None,
                      variants_only: bool, extra_mpileup: str, extra_call: str) -> str:
    pile = ["bcftools mpileup", f"-f {q(ref)}", f"-a {q(annotations)}"]
    if regions:
        pile.append(f"-R {q(regions)}")
    if ignore_rg:
        pile.append("--ignore-RG")
    if skip_indels:
        pile.append("-I")
    if max_depth is not None:
        pile.append(f"-d {int(max_depth)}")
    if min_mapq is not None:
        pile.append(f"-q {int(min_mapq)}")
    if min_baseq is not None:
        pile.append(f"-Q {int(min_baseq)}")
    if extra_mpileup:
        pile.append(extra_mpileup)
    if bam_list:
        pile.append(f"--bam-list {q(bam_list)}")
    else:
        pile.extend(q(b) for b in (bams or []))
    pile.append("-Ou")

    call = ["bcftools call", "-m", f"--ploidy {q(ploidy)}"]
    if variants_only:
        call.append("-v")
    if extra_call:
        call.append(extra_call)
    call.append(f"-O{out_flag(out)} -o {q(out)}")
    return " ".join(pile) + " | " + " ".join(call)


#: What ``--bam-dir`` picks up. CRAM as well as BAM, because ``--bam`` takes either and a
#: directory of CRAMs reporting "no alignments" would be a silly way to find that out.
BAM_DIR_PATTERNS = ("*.bam", "*.cram")


def bams_in_dir(directory: str, patterns: tuple[str, ...] = BAM_DIR_PATTERNS) -> list[str]:
    """Every alignment directly in ``directory``, sorted.

    Sorted because the order decides the sample order of the output callset, and glob order
    is whatever the filesystem says -- two runs over the same directory should not produce
    columns in different orders. Not recursive: a directory of alignments is what this is
    for, and walking a tree would quietly pick up an unrelated subdirectory of BAMs.
    """
    if not os.path.isdir(directory):
        raise SystemExit(f"call_variants: --bam-dir {directory} is not a directory")
    found = sorted(str(p) for pat in patterns
                   for p in pathlib.Path(directory).glob(pat))
    if not found:
        raise SystemExit(f"call_variants: no {' / '.join(patterns)} in {directory}")
    return found


def call_variants(ref: str, out: str, *, bams: list[str] | None = None,
                  bam_list: str | None = None, bam_dir: str | None = None,
                  regions: str | None = None,
                  threads: int = 1, chunk_size: int | None = None,
                  annotations: str = BCFTOOLS_MPILEUP_ANNOTATIONS, ploidy: str = "2",
                  ignore_rg: bool = True, sample_suffix: str | None = ".bam",
                  skip_indels: bool = False,
                  max_depth: int | None = None, min_mapq: int | None = None,
                  min_baseq: int | None = None, variants_only: bool = False,
                  extra_mpileup: str = "",
                  extra_call: str = "", keep_chunks: str | None = None,
                  dry_run: bool = False) -> list[str]:
    """Call variants, splitting a region list over ``threads`` concurrent jobs.

    Parameters
    ----------
    bams, bam_list, bam_dir:
        The alignments, given one of three ways: ``bams`` as paths, ``bam_list`` as a file
        of paths one per line (what `bcftools mpileup --bam-list` reads), or ``bam_dir`` as
        a directory whose alignments are used, sorted -- see :func:`bams_in_dir`.
    regions:
        Optional region file. A ``.bed`` is read 0-based half-open, anything else as
        1-based ``CHROM POS``, which is bcftools' own rule and is preserved when splitting.
    threads:
        Concurrent chunk jobs. This is process-level parallelism over the region list --
        ``bcftools mpileup --threads`` only parallelises compression, so it is not what
        makes calling faster.
    chunk_size:
        Regions per chunk. The default splits the list into ``threads`` pieces, which is
        the fewest jobs that still keeps every thread busy; set it to make chunks a fixed
        size instead (smaller chunks even out uneven regions at the cost of more jobs).
    ignore_rg, sample_suffix:
        One BAM per sample is the usual arrangement, so ``--ignore-RG`` is on by default:
        each alignment becomes one sample regardless of what its read groups say. bcftools
        then names each sample after the *path* it was given, so the samples are renamed
        afterwards to the file name with ``sample_suffix`` removed (``.bam`` by default).
        Pass the whole trailing part to strip -- ``.sorted.dup.pf.bam`` -- or ``None`` to
        keep the names bcftools wrote. A suffix that would make two samples collide is an
        error, raised before any calling starts rather than after.
    variants_only:
        Emit only variant sites (`bcftools call -v`). Off by default, because calling a
        list of known positions is usually about filling them in: a reference call at a
        target position is the answer "this sample is reference here", and `-v` would drop
        it. Turn it on for whole-genome calling, where the non-variant sites are just bulk.
    keep_chunks:
        Directory to leave the per-chunk BCFs and region files in. They go to a temporary
        directory otherwise, and are removed once concatenated.
    dry_run:
        Return the commands without running any of them.

    Returns the commands run, in order, so a run can be reproduced or inspected.
    """
    given = [n for n, v in (("--bam", bams), ("--bam-list", bam_list),
                            ("--bam-dir", bam_dir)) if v]
    if not given:
        raise SystemExit("call_variants: give --bam (one or more), --bam-list or --bam-dir")
    if len(given) > 1:
        raise SystemExit(f"call_variants: {', '.join(given)} are alternatives, not both")
    if bam_dir:
        bams = bams_in_dir(bam_dir)
        print(f"  {len(bams)} alignment(s) in {bam_dir}")
    if threads < 1:
        raise SystemExit("call_variants: --threads must be at least 1")
    if not dry_run:
        require("bcftools")

    # fail on a naming clash before spending an hour calling, not after
    if ignore_rg and sample_suffix is not None:
        paths = bams
        if paths is None:
            if os.path.exists(bam_list):
                paths = [ln.strip() for ln in open(bam_list) if ln.strip()]
            elif not dry_run:
                raise SystemExit(f"call_variants: --bam-list {bam_list} does not exist")
        if paths:
            sample_rename_map(paths, sample_suffix)

    common = dict(ref=ref, bams=bams, bam_list=bam_list, annotations=annotations,
                  ploidy=ploidy, ignore_rg=ignore_rg, skip_indels=skip_indels,
                  max_depth=max_depth, min_mapq=min_mapq, min_baseq=min_baseq,
                  variants_only=variants_only, extra_mpileup=extra_mpileup,
                  extra_call=extra_call)

    # nothing to split: one region file-less job, or one chunk's worth of regions
    workdir = keep_chunks or tempfile.mkdtemp(prefix="call_variants.")
    if keep_chunks:
        os.makedirs(keep_chunks, exist_ok=True)
    try:
        chunks = []
        if regions and threads > 1:
            chunks = split_regions(regions, workdir, n_chunks=threads,
                                   chunk_size=chunk_size)
        if not chunks:
            cmd = _mpileup_call_cmd(out=out, regions=regions, **common)
            if dry_run:
                return [cmd] + _rename_cmds(ignore_rg, sample_suffix, out)
            sh(cmd, tools=("bcftools",))
            index_vcf(out)
            if ignore_rg and sample_suffix is not None:
                _rename_samples(out, sample_suffix)
            return [cmd]

        parts = [os.path.join(workdir, f"part{i:04d}.bcf") for i in range(len(chunks))]
        cmds = [_mpileup_call_cmd(out=p, regions=c, **common)
                for p, c in zip(parts, chunks)]
        # each part is indexed so concat can merge them by coordinate rather than trusting
        # the order the region file happened to be in
        cmds += [f"bcftools index {q(p)}" for p in parts]
        concat = (f"bcftools concat -a {' '.join(q(p) for p in parts)} "
                  f"-O{out_flag(out)} -o {q(out)}")
        cmds.append(concat)
        if dry_run:
            return cmds + _rename_cmds(ignore_rg, sample_suffix, out)

        with ThreadPoolExecutor(max_workers=threads) as pool:
            list(pool.map(lambda c: sh(c, tools=("bcftools",)), cmds[:len(parts)]))
        for p in parts:
            index_vcf(p)
        sh(concat, tools=("bcftools",))
        index_vcf(out)
        if ignore_rg and sample_suffix is not None:
            _rename_samples(out, sample_suffix)
        return cmds
    finally:
        if not keep_chunks and os.path.isdir(workdir):
            shutil.rmtree(workdir, ignore_errors=True)


def n_regions(path: str) -> int:
    """How many regions a file lists, for reporting how the work was split."""
    return len(_region_lines(path))
