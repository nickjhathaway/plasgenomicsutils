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

**Hundreds of BAMs on a network filesystem need splitting the other way too.** Every
region job opens every alignment, so 600 BAMs over 5 jobs is 6,000 concurrent file
handles (BAM plus index) against one NFS server, and the run stalls in I/O rather than
computing. ``bam_batch`` caps that the way the CNV pipeline's Fws step does: the
alignments are cut into contiguous groups of at most that many, one job runs per (group x
region chunk), and the peak is ``threads x bam_batch`` handles whatever the cohort size.
The price is that the groups are called separately: a group where no sample carries an
ALT does not emit it, so the group callsets disagree on alleles. :mod:`.harmonize` puts
them back on one allele set -- an allele a group never saw gets AD 0 for its samples,
which is what that group's pileup found -- and ``bcftools merge`` joins them. Genotypes,
AD, ADF, ADR and the per-site read counts (INFO/AD, ADF, ADR, DP, DP4) come out exactly
as the joint call would have them; what cannot be recovered are the per-site *statistics*
bcftools computes over the pooled reads -- the ``*BZ`` z-scores, ``FS``, ``MQ``, QUAL --
which are combined across groups by rule instead (:data:`GROUP_MERGE_RULES`). Harmonizing
is SNP-only, so indel records are dropped in this mode, and PL goes with them.
"""

from __future__ import annotations

import math
import os
import pathlib
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

from . import harmonize as H
from .bcftools import index_vcf, info_tags, out_flag, q, require, sample_names, sh
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


def bam_groups(paths: list[str], bam_batch: int | None) -> list[list[str]]:
    """Cut the alignments into contiguous groups of at most ``bam_batch``.

    Contiguous, so the merged callset keeps the input sample order (``bcftools merge``
    appends the files' samples in order). One group -- ``bam_batch`` off, or at least as
    large as the cohort -- means nothing to harmonize, and calling runs as it would have
    without it.
    """
    if not bam_batch or bam_batch >= len(paths):
        return [list(paths)]
    size = int(bam_batch)
    return [list(paths[i:i + size]) for i in range(0, len(paths), size)]


#: How ``bcftools merge -i`` combines each INFO tag across the per-group callsets. The
#: counts are sums, which is exactly what the joint call would have counted. The
#: statistics are a stand-in: ``FS`` is a p-value, so the smallest (most significant) is
#: kept; the Mann-Whitney z-scores and the quality summaries are averaged. A z pooled over
#: all the reads would grow with the number of groups (~sqrt(G) for the same effect), so
#: these are conservative in the *keep* direction for ``hard_qc_filter`` -- its effect-size
#: guard (``bias_eff``, from ADF/ADR) is on the exact counts and unaffected. Any tag not in
#: the header is left out of the rule at run time; the merge default for the rest is the
#: first group's value.
GROUP_MERGE_RULES = {
    "DP": "sum", "DP4": "sum", "AD": "sum", "ADF": "sum", "ADR": "sum", "SCR": "sum",
    "FS": "min",
    "RPBZ": "avg", "MQBZ": "avg", "MQSBZ": "avg", "BQBZ": "avg", "SCBZ": "avg",
    "MQ": "avg", "MQ0F": "avg", "VDB": "avg", "SGB": "avg",
}

#: INFO fields dropped before the groups are merged: ``Number=A`` and count-of-alleles
#: fields ``bcftools merge`` recomputes from the merged genotypes anyway.
GROUP_STALE_INFO = ("AC", "AN", "AF")

#: FORMAT fields kept through harmonizing in group mode. The ``Number=R`` counts are
#: re-laid-out on the union allele set with zeros; everything else per-allele (PL, which is
#: ``Number=G``) is dropped, since a likelihood for a genotype the group never scored does
#: not exist.
GROUP_KEEP_FORMAT = ("GT", "AD", "ADF", "ADR")


def merge_info_rules(path: str, rules: dict = GROUP_MERGE_RULES) -> str:
    """The ``bcftools merge -i`` argument for a file, restricted to tags its header has."""
    present = info_tags(path)
    return ",".join(f"{t}:{r}" for t, r in rules.items() if t in present)


def _merge_cmd(harmonized: list[str], out: str, rules: str) -> str:
    return (f"bcftools merge -i {q(rules)} {' '.join(q(p) for p in harmonized)} "
            f"-O{out_flag(out)} -o {q(out)}")


def _harmonize_groups(group_files: list[str], workdir: str, *,
                      keep_ref_only: bool) -> tuple[list[str], list[str]]:
    """Put the per-group callsets on one allele set so ``bcftools merge`` can join them.

    No cleaning (``min_ad=0``, ``min_af=0``): the alleles are what bcftools called, and
    the only edit is adding, with zero depth, the ones a group did not see. Genotypes are
    re-indexed, not re-called. Returns the harmonized BCF paths and the commands run.
    """
    union, _dups, ambiguous, st = H.accumulate_union(
        group_files, 0, 0.0, 0.2, drop_indels=True, keep_ref_only=keep_ref_only)
    if ambiguous:
        print(f"  WARNING: {len(ambiguous)} position(s) carried more than one record with "
              f"ALTs in a group callset; the one with most ALTs was kept")
    stale_fmt = sorted(set().union(
        *(set(H.stale_format_fields(f, keep=GROUP_KEEP_FORMAT)) for f in group_files)))
    strip = ("-x " + q(",".join("FORMAT/" + f for f in stale_fmt)) + " ") if stale_fmt else ""
    n_indel = sum(v["indel_context"] for v in st["per_file"].values())
    print(f"  harmonizing {len(group_files)} group callsets: {st['union_sites']} site(s), "
          f"{st['union_with_alts']} with an ALT in some group"
          + (f", {n_indel} indel-context record(s) dropped" if n_indel else "")
          + (f"; dropping FORMAT/{','.join(stale_fmt)}" if stale_fmt else ""))
    out_files, cmds, absent = [], [], 0
    for i, f in enumerate(group_files):
        tmp = os.path.join(workdir, f"h_group{i:04d}.tmp.vcf")
        h = os.path.join(workdir, f"h_group{i:04d}.bcf")
        r = H.harmonize_file(f, tmp, union, 0, 0.0, 0.2, drop_indels=True,
                             regenotype=False, stale_info=GROUP_STALE_INFO)
        absent = max(absent, r["absent"])
        cmd = f"bcftools annotate {strip}-Ob -o {q(h)} {q(tmp)}"
        sh(cmd, tools=("bcftools",))
        cmds.append(cmd)
        os.remove(tmp)
        index_vcf(h)
        out_files.append(h)
    if absent:
        print(f"  NOTE: up to {absent} site(s) were emitted by some groups but not others "
              f"(--variants-only), so those samples get missing genotypes after the merge")
    return out_files, cmds


def _harmonize_dry_run_cmds(group_files: list[str], workdir: str) -> list[str]:
    hs = [os.path.join(workdir, f"h_group{i:04d}.bcf") for i in range(len(group_files))]
    lines = [f"# harmonize (in-process): union the ALTs of {' '.join(q(g) for g in group_files)}, "
             f"zero-fill AD/ADF/ADR for alleles a group did not see, re-index GT, "
             f"drop indel records and INFO/{','.join(GROUP_STALE_INFO)}"]
    lines += [f"bcftools annotate -x <per-genotype FORMAT fields, e.g. FORMAT/PL> -Ob "
              f"-o {q(h)} {q(h[:-4] + '.tmp.vcf')}" for h in hs]
    return lines


def _read_bam_list(bam_list: str) -> list[str]:
    with open(bam_list) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


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
                  bam_batch: int | None = None,
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
    bam_batch:
        Alignments per group. Off (``None``/0) calls every alignment in every job. Set it
        -- 100 is what the CNV pipeline uses -- and the alignments are cut into contiguous
        groups of at most that many, one job runs per group and region chunk, and the
        group callsets are harmonized (:mod:`.harmonize`, no cleaning, genotypes
        re-indexed rather than re-called) and ``bcftools merge``-d with
        :data:`GROUP_MERGE_RULES`. That caps concurrently open alignments at
        ``threads x bam_batch``, which is what keeps a large cohort on NFS from stalling.
        A cohort no larger than the batch is one group and calls exactly as without it.
        See the module docstring for what the split changes: SNP-only output, no PL, and
        per-site bias statistics combined by rule rather than computed over all reads.
        Because it is SNP-only, ``skip_indels`` is required whenever the split is in
        effect -- an error rather than a note, so indel records never vanish unannounced.
    keep_chunks:
        Directory to leave the per-chunk BCFs and region files in -- and in group mode the
        per-group callsets, harmonized and not. They go to a temporary directory
        otherwise, and are removed once concatenated.
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

    paths = bams
    if paths is None:
        if os.path.exists(bam_list):
            paths = _read_bam_list(bam_list)
        elif not dry_run:
            raise SystemExit(f"call_variants: --bam-list {bam_list} does not exist")
    # fail on a naming clash before spending an hour calling, not after
    if ignore_rg and sample_suffix is not None and paths:
        sample_rename_map(paths, sample_suffix)

    if bam_batch is not None and bam_batch < 0:
        raise SystemExit("call_variants: --bam-batch must be 0 (off) or a positive number")
    if bam_batch and paths is None:
        raise SystemExit(f"call_variants: --bam-batch needs the alignments listed, and "
                         f"--bam-list {bam_list} does not exist")
    if bam_batch and paths and len(paths) > bam_batch and not skip_indels:
        # Harmonizing the groups is SNP-only, so indel records would be called and then
        # silently dropped. Make the caller say so, so nobody wonders where they went.
        raise SystemExit(
            "call_variants: --bam-batch output is SNP-only (indel records are dropped when "
            "the group callsets are harmonized), so it requires --skip-indels: pass it to "
            "confirm, and the groups do not spend time computing indels either")
    groups = bam_groups(paths, bam_batch) if paths else [None]
    if len(groups) == 1 and paths and not bam_batch and len(paths) * threads > 1000:
        print(f"  NOTE: {len(paths)} alignments x {threads} jobs = "
              f"{len(paths) * threads} concurrently open alignments (plus indexes). On a "
              f"network filesystem that can stall; --bam-batch 100 caps it at "
              f"{100 * threads}. See --help.")

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
        if len(groups) > 1:
            return _call_in_groups(
                groups, chunks, regions=regions, out=out, workdir=workdir,
                threads=threads, common=common, ignore_rg=ignore_rg,
                sample_suffix=sample_suffix, variants_only=variants_only, dry_run=dry_run)
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


def _call_in_groups(groups: list[list[str]], chunks: list[str], *, regions: str | None,
                    out: str, workdir: str, threads: int, common: dict, ignore_rg: bool,
                    sample_suffix: str | None, variants_only: bool,
                    dry_run: bool) -> list[str]:
    """Group mode: a job per (alignment group x region chunk), then harmonize and merge.

    Jobs are ordered group-major, so the ``threads`` running at any moment are mostly
    over the same ``bam_batch`` files rather than ``threads`` different sets of them --
    the gentlest pattern for a network filesystem's cache.
    """
    region_files = chunks or [regions]
    n_bams = sum(len(g) for g in groups)
    print(f"  {n_bams} alignment(s) in {len(groups)} group(s) of <= {max(map(len, groups))} "
          f"x {len(region_files)} region chunk(s) = {len(groups) * len(region_files)} "
          f"job(s) over {threads} thread(s); at most ~{threads * max(map(len, groups))} "
          f"alignments open at once")
    print("  group mode: SNP-only output, no PL; per-site bias statistics are combined "
          "across groups by rule, not computed over all reads -- see call_variants --help")

    common = {k: v for k, v in common.items() if k not in ('bams', 'bam_list')}
    list_files, group_files, tiles, cmds = [], [], [], []
    for gi, group in enumerate(groups):
        lst = os.path.join(workdir, f"group{gi:04d}.bams.txt")
        list_files.append(lst)
        if not dry_run:
            with open(lst, "w") as fh:
                fh.writelines(p + "\n" for p in group)
        gout = os.path.join(workdir, f"group{gi:04d}.bcf")
        group_files.append(gout)
        if len(region_files) == 1:
            tiles.append((gi, None, gout))
        else:
            for ci in range(len(region_files)):
                tiles.append((gi, ci, os.path.join(workdir, f"group{gi:04d}.part{ci:04d}.bcf")))
    tile_cmds = [_mpileup_call_cmd(out=p, regions=region_files[ci or 0], bams=None,
                                   bam_list=list_files[gi], **common)
                 for gi, ci, p in tiles]
    cmds += tile_cmds
    concat_cmds = []
    if len(region_files) > 1:
        for gi, gout in enumerate(group_files):
            parts = [p for g, _c, p in tiles if g == gi]
            cmds += [f"bcftools index {q(p)}" for p in parts]
            concat_cmds.append(f"bcftools concat -a {' '.join(q(p) for p in parts)} "
                               f"-Ob -o {q(gout)}")
        cmds += concat_cmds
    harmonized = [os.path.join(workdir, f"h_group{gi:04d}.bcf") for gi in range(len(groups))]
    if dry_run:
        rules = ",".join(f"{t}:{r}" for t, r in GROUP_MERGE_RULES.items())
        return (cmds + _harmonize_dry_run_cmds(group_files, workdir)
                + [_merge_cmd(harmonized, out, rules)]
                + _rename_cmds(ignore_rg, sample_suffix, out))

    with ThreadPoolExecutor(max_workers=threads) as pool:
        list(pool.map(lambda c: sh(c, tools=("bcftools",)), tile_cmds))
    if len(region_files) > 1:
        for _gi, _ci, p in tiles:
            index_vcf(p)
        with ThreadPoolExecutor(max_workers=threads) as pool:
            list(pool.map(lambda c: sh(c, tools=("bcftools",)), concat_cmds))
    for gout in group_files:
        index_vcf(gout)

    harmonized, hcmds = _harmonize_groups(group_files, workdir,
                                          keep_ref_only=not variants_only)
    cmds += hcmds
    merge = _merge_cmd(harmonized, out, merge_info_rules(harmonized[0]))
    sh(merge, tools=("bcftools",))
    cmds.append(merge)
    index_vcf(out)
    if ignore_rg and sample_suffix is not None:
        _rename_samples(out, sample_suffix)
    return cmds


def n_regions(path: str) -> int:
    """How many regions a file lists, for reporting how the work was split."""
    return len(_region_lines(path))
