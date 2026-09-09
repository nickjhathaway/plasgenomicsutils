"""call_variants: the annotations match the filter, and splitting regions changes nothing.

Calling is one core per job, so the speed-up comes from splitting the region list rather
than from bcftools' own --threads. What has to hold is that the split is invisible in the
result, and that what gets called carries what `hard_qc_filter --caller bcftools` reads.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from plasgenomicsutils.lib import call_variants as C
from plasgenomicsutils.lib.vcf_filters import (BCFTOOLS_MPILEUP_ANNOTATIONS,
                                               BCFTOOLS_MPILEUP_DEFAULT_TAGS,
                                               _BCFTOOLS_QC_NEEDS)

needs_bcftools = pytest.mark.skipif(not shutil.which("bcftools"),
                                    reason="bcftools not on PATH")


def _bed(tmp_path, n, name="sites.bed"):
    p = tmp_path / name
    p.write_text("".join(f"chr1\t{1000 + 100 * i}\t{1001 + 100 * i}\n" for i in range(n)))
    return str(p)


# ---- the two lists cannot drift apart --------------------------------------------

def test_everything_the_qc_filter_reads_is_obtainable_from_the_default_call():
    """The guard this wrapper exists for.

    A bcftools comparison against an absent tag is false, so a callset made without the
    right -a list passes hard_qc_filter untouched. Calling through here has to request
    every INFO tag the filter tests that mpileup does not already write.
    """
    requested = {a.split("/")[-1] for a in BCFTOOLS_MPILEUP_ANNOTATIONS.split(",")
                 if a.startswith("INFO/")}
    needed = {t for tags in _BCFTOOLS_QC_NEEDS.values() for t in tags}
    missing = needed - requested - set(BCFTOOLS_MPILEUP_DEFAULT_TAGS)
    assert not missing, (
        f"hard_qc_filter --caller bcftools reads INFO/{', INFO/'.join(sorted(missing))}, "
        "which call_variants neither requests nor gets for free")


# ---- splitting a region file -----------------------------------------------------

def test_chunks_keep_the_extension_so_the_coordinates_keep_their_meaning(tmp_path):
    # bcftools reads 0-based BED vs 1-based CHROM/POS off the extension, so a chunk of a
    # .bed that is not a .bed silently shifts every region by one
    out = tmp_path / "chunks"
    out.mkdir()
    chunks = C.split_regions(_bed(tmp_path, 8), str(out), n_chunks=4)
    assert len(chunks) == 4
    assert all(c.endswith(".bed") for c in chunks)
    assert sum(len(open(c).readlines()) for c in chunks) == 8


def test_a_double_extension_survives_too(tmp_path):
    src = tmp_path / "sites.bed.gz"          # name only; contents are plain here
    src.write_text("chr1\t1\t2\nchr1\t3\t4\n")
    out = tmp_path / "c"
    out.mkdir()
    chunks = C.split_regions(str(src), str(out), n_chunks=2)
    assert all(c.endswith(".bed.gz") for c in chunks)


def test_threads_decides_the_number_of_chunks_and_chunk_size_overrides_it(tmp_path):
    out = tmp_path / "c"
    out.mkdir()
    # 400 regions over 8 jobs is 50 apiece -- the split the recipe this replaces did by hand
    assert len(C.split_regions(_bed(tmp_path, 400), str(out), n_chunks=8)) == 8
    assert len(open(C.split_regions(_bed(tmp_path, 400), str(out), n_chunks=8)[0])
               .readlines()) == 50
    # an explicit size wins
    chunks = C.split_regions(_bed(tmp_path, 400), str(out), n_chunks=8, chunk_size=25)
    assert len(chunks) == 16


def test_one_region_is_not_split(tmp_path):
    out = tmp_path / "c"
    out.mkdir()
    assert C.split_regions(_bed(tmp_path, 1), str(out), n_chunks=8) == []


def test_comments_and_blank_lines_are_not_counted_as_regions(tmp_path):
    p = tmp_path / "sites.bed"
    p.write_text("# a header\ntrack name=x\n\nchr1\t1\t2\nchr1\t3\t4\n")
    assert C.n_regions(str(p)) == 2


def test_an_empty_region_file_is_an_error(tmp_path):
    p = tmp_path / "empty.bed"
    p.write_text("# nothing here\n")
    with pytest.raises(SystemExit, match="lists no regions"):
        C.split_regions(str(p), str(tmp_path), n_chunks=2)


# ---- the commands built ----------------------------------------------------------

def test_dry_run_builds_the_split_call_index_concat_sequence(tmp_path):
    cmds = C.call_variants("ref.fa", str(tmp_path / "out.bcf"), bams=["a.bam"],
                           regions=_bed(tmp_path, 6), threads=3, dry_run=True)
    calls = [c for c in cmds if "mpileup" in c]
    assert len(calls) == 3                                  # one job per chunk
    assert all("-R " in c and ".bed" in c for c in calls)
    assert all(f"-a {BCFTOOLS_MPILEUP_ANNOTATIONS}" in c for c in calls)
    assert sum("bcftools index" in c for c in cmds) == 3     # concat needs them indexed
    assert any(c.startswith("bcftools concat -a ") for c in cmds)
    assert "reheader" in cmds[-1]                            # naming is the last step


def test_without_regions_it_is_a_single_job(tmp_path):
    cmds = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bams=["a.bam"], threads=8,
                           dry_run=True)
    calls = [c for c in cmds if "mpileup" in c]
    assert len(calls) == 1 and "-R " not in calls[0]         # nothing to split
    assert not any("concat" in c for c in cmds)


def test_bam_list_and_the_passthrough_options(tmp_path):
    cmd = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bam_list="bams.txt",
                          ploidy="1", ignore_rg=True, skip_indels=True, max_depth=500,
                          min_mapq=20, extra_mpileup="--no-BAQ",
                          extra_call="--variants-only", dry_run=True)[0]
    assert "--bam-list bams.txt" in cmd and "--ignore-RG" in cmd
    assert " -I " in cmd and "-d 500" in cmd and "-q 20" in cmd and "--no-BAQ" in cmd
    assert "--ploidy 1" in cmd and "--variants-only" in cmd


def test_the_alignment_arguments_are_required_and_exclusive(tmp_path):
    with pytest.raises(SystemExit, match="--bam-list or --bam-dir"):
        C.call_variants("ref.fa", "o.bcf", dry_run=True)
    for kw in ({"bams": ["a.bam"], "bam_list": "b.txt"},
               {"bams": ["a.bam"], "bam_dir": "d"},
               {"bam_list": "b.txt", "bam_dir": "d"}):
        with pytest.raises(SystemExit, match="alternatives, not both"):
            C.call_variants("ref.fa", "o.bcf", dry_run=True, **kw)


# ---- naming samples after their files --------------------------------------------

def test_the_default_strips_bam_from_the_file_name(tmp_path):
    # mpileup --ignore-RG names each sample after the path it was handed
    assert C.sample_name_for("/tank/wgs/17017-227227.bam") == "17017-227227"
    assert C.sample_name_for("/tank/wgs/17017-227227.sorted.dup.pf.bam") == \
        "17017-227227.sorted.dup.pf"
    # the whole pipeline suffix, which is what a real BAM tends to carry
    assert C.sample_name_for("/tank/wgs/17017-227227.sorted.dup.pf.bam",
                             ".sorted.dup.pf.bam") == "17017-227227"
    # a suffix that is not there leaves the name alone
    assert C.sample_name_for("/tank/wgs/s1.cram") == "s1.cram"
    assert C.sample_name_for("/tank/wgs/s1.bam", None) == "s1.bam"


def test_a_suffix_that_would_collide_is_refused_before_anything_is_called(tmp_path):
    # same file name under two run directories: renaming would merge two samples into one
    with pytest.raises(SystemExit, match="same sample name"):
        C.sample_rename_map(["runA/s1.bam", "runB/s1.bam"])
    # and the check runs up front, so an hour of calling is not spent first
    with pytest.raises(SystemExit, match="same sample name"):
        C.call_variants("ref.fa", str(tmp_path / "o.bcf"),
                        bams=["runA/s1.bam", "runB/s1.bam"], dry_run=True)
    # keeping them apart is enough
    assert C.sample_rename_map(["runA/s1.bam", "runB/s2.bam"]) == {
        "runA/s1.bam": "s1", "runB/s2.bam": "s2"}


def test_a_suffix_that_consumes_the_whole_name_is_refused():
    with pytest.raises(SystemExit, match="leaves nothing to call the sample"):
        C.sample_rename_map(["/x/s1.bam"], "s1.bam")


def test_ignore_rg_is_on_by_default_and_can_be_turned_off(tmp_path):
    on = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bams=["a.bam"], dry_run=True)
    assert "--ignore-RG" in on[0]
    assert any("reheader" in c for c in on)          # and the rename follows
    off = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bams=["a.bam"],
                          ignore_rg=False, dry_run=True)
    assert "--ignore-RG" not in off[0]
    assert not any("reheader" in c for c in off)     # RG names are already real names
    keep = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bams=["a.bam"],
                           sample_suffix=None, dry_run=True)
    assert not any("reheader" in c for c in keep)


# ---- end to end ------------------------------------------------------------------

@needs_bcftools
def test_splitting_the_regions_does_not_change_the_calls(tmp_path):
    """The property the parallelism has to have.

    Same positions, same genotypes, same depths and bias statistics. QUAL is left out of
    the comparison on purpose: `bcftools mpileup -R` computes indel likelihoods and BAQ
    from the reads around a position, so which neighbours share a chunk moves QUAL by a
    few points at indels and the occasional SNP next to one. That is bcftools' own
    behaviour -- cutting a region file in half by hand reproduces it -- and nothing in
    this package reads QUAL.
    """
    pysam = pytest.importorskip("pysam")
    import random
    random.seed(3)
    ref = "".join(random.choice("ACGT") for _ in range(3000))
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
    pysam.faidx(str(fa))

    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"LN": len(ref), "SN": "chr1"}],
           "RG": [{"ID": "rg1", "SM": "s1", "PL": "ILLUMINA"}]}
    bam = tmp_path / "in.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        for site in (1000, 1500, 2000):
            for i in range(30):
                st = site - 50
                seq = list(ref[st:st + 100])
                if i % 2 == 0:
                    seq[50] = "T" if ref[site] != "T" else "G"
                a = pysam.AlignedSegment()
                a.query_name = f"r{site}_{i}"; a.query_sequence = "".join(seq)
                a.flag = 16 if i % 4 in (2, 3) else 0
                a.reference_id = 0; a.reference_start = st; a.mapping_quality = 60
                a.cigarstring = "100M"
                a.query_qualities = pysam.qualitystring_to_array("I" * 100)
                a.set_tag("RG", "rg1")
                out.write(a)
    pysam.index(str(bam))

    beds = tmp_path / "sites.bed"
    beds.write_text("chr1\t1000\t1001\nchr1\t1500\t1501\nchr1\t2000\t2001\n")

    serial = str(tmp_path / "serial.bcf")
    par = str(tmp_path / "par.bcf")
    C.call_variants(str(fa), serial, bams=[str(bam)], regions=str(beds), threads=1)
    C.call_variants(str(fa), par, bams=[str(bam)], regions=str(beds), threads=3)

    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%DP\t%FS\t%RPBZ\t%MQBZ[\t%AD\t%GT]\n"
    rows = [subprocess.run(["bcftools", "query", "-f", fmt, p], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, text=True).stdout for p in (serial, par)]
    assert rows[0] == rows[1]
    assert rows[0].count("\n") == 3


@needs_bcftools
def test_samples_come_out_named_after_their_files(tmp_path):
    """The default path: --ignore-RG plus the rename, end to end."""
    pysam = pytest.importorskip("pysam")
    import random
    from plasgenomicsutils.lib.bcftools import sample_names as names_of
    random.seed(5)
    ref = "".join(random.choice("ACGT") for _ in range(1500))
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
    pysam.faidx(str(fa))

    bams = []
    for sample in ("17017-227227", "90909-962962"):
        # read groups say something else entirely; --ignore-RG means the file wins
        hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
               "SQ": [{"LN": len(ref), "SN": "chr1"}],
               "RG": [{"ID": "rg1", "SM": "WRONG_NAME", "PL": "ILLUMINA"}]}
        path = tmp_path / f"{sample}.sorted.dup.pf.bam"
        with pysam.AlignmentFile(str(path), "wb", header=hdr) as out:
            for i in range(20):
                seq = list(ref[700:800])
                if i % 2 == 0:
                    seq[50] = "T" if ref[750] != "T" else "G"
                a = pysam.AlignedSegment()
                a.query_name = f"r{i}"; a.query_sequence = "".join(seq)
                a.flag = 16 if i % 4 in (2, 3) else 0
                a.reference_id = 0; a.reference_start = 700; a.mapping_quality = 60
                a.cigarstring = "100M"
                a.query_qualities = pysam.qualitystring_to_array("I" * 100)
                a.set_tag("RG", "rg1")
                out.write(a)
        pysam.index(str(path))
        bams.append(str(path))

    out = str(tmp_path / "out.bcf")
    C.call_variants(str(fa), out, bams=bams)
    assert names_of(out) == ["17017-227227.sorted.dup.pf", "90909-962962.sorted.dup.pf"]

    out2 = str(tmp_path / "out2.bcf")
    C.call_variants(str(fa), out2, bams=bams, sample_suffix=".sorted.dup.pf.bam")
    assert names_of(out2) == ["17017-227227", "90909-962962"]

    # the RG tags said WRONG_NAME; --no-ignore-rg is what would have used them
    out3 = str(tmp_path / "out3.bcf")
    C.call_variants(str(fa), out3, bams=bams, ignore_rg=False)
    assert names_of(out3) == ["WRONG_NAME"]


def test_the_default_pipeline_config_makes_the_caller_switchable():
    """`--emit-default-config` writes `caller` out, like `keep_bed`.

    A bcftools callset carries none of GATK's metrics, so which vocabulary the step reads
    is the setting someone needs to find; JSON has no comments, so a key that is not
    written is a key nobody knows exists.
    """
    from plasgenomicsutils.lib.filter_pipeline import DEFAULT_CONFIG
    step = next(s for s in DEFAULT_CONFIG["steps"] if s["name"] == "hard_qc_filter")
    assert step.get("params", {}).get("caller") == "gatk"


# ---- --bam-dir: a directory of alignments instead of a list ------------------------

def _bam_dir(tmp_path, names=("b.bam", "a.bam", "c.cram", "notes.txt", "a.bam.bai")):
    d = tmp_path / "bams"
    d.mkdir()
    for n in names:
        (d / n).write_text("")
    (d / "nested").mkdir()
    (d / "nested" / "deep.bam").write_text("")
    return str(d)


def test_bams_in_dir_is_sorted_and_not_recursive(tmp_path):
    from plasgenomicsutils.lib.call_variants import bams_in_dir

    got = [os.path.basename(p) for p in bams_in_dir(_bam_dir(tmp_path))]
    # sorted, because this order becomes the sample order of the output callset
    assert got == ["a.bam", "b.bam", "c.cram"]
    # a subdirectory of alignments is somebody else's cohort, not part of this one
    assert not any("deep" in p for p in got)
    # and nothing that is not an alignment
    assert not any(p.endswith((".txt", ".bai")) for p in got)


def test_an_empty_or_missing_directory_says_which(tmp_path):
    from plasgenomicsutils.lib.call_variants import bams_in_dir

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(SystemExit, match="no \\*.bam"):
        bams_in_dir(str(empty))
    with pytest.raises(SystemExit, match="is not a directory"):
        bams_in_dir(str(tmp_path / "nope"))


def test_bam_dir_gives_the_same_commands_as_listing_those_files(tmp_path):
    d = _bam_dir(tmp_path)
    listed = sorted(str(p) for p in Path(d).glob("*.bam")) + [str(Path(d) / "c.cram")]
    a = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bam_dir=d, dry_run=True)
    b = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bams=listed, dry_run=True)
    assert a == b


def test_a_directory_still_gets_its_samples_renamed(tmp_path):
    """The rename runs however the alignments arrived -- the directory form must not slip
    past it and leave samples named after their paths."""
    d = tmp_path / "cohort"
    d.mkdir()
    for n in ("s1.sorted.bam", "s2.sorted.bam"):
        (d / n).write_text("")
    cmds = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bam_dir=str(d),
                           sample_suffix=".sorted.bam", dry_run=True)
    assert any("bcftools reheader" in c and ".sorted.bam" in c for c in cmds)


# ---- splitting the alignments, for cohorts too big to open all at once ------------

def test_bam_groups_are_contiguous_and_off_below_the_batch():
    paths = [f"s{i}.bam" for i in range(7)]
    assert C.bam_groups(paths, None) == [paths]
    assert C.bam_groups(paths, 0) == [paths]
    assert C.bam_groups(paths, 7) == [paths]                 # one group: nothing to merge
    assert C.bam_groups(paths, 3) == [paths[:3], paths[3:6], paths[6:]]


def test_group_mode_dry_run_tiles_groups_by_chunks_then_harmonizes_and_merges(tmp_path):
    bams = [f"s{i}.bam" for i in range(5)]
    cmds = C.call_variants("ref.fa", str(tmp_path / "out.bcf"), bams=bams,
                           regions=_bed(tmp_path, 6), threads=3, bam_batch=2, skip_indels=True,
                           dry_run=True)
    calls = [c for c in cmds if "mpileup" in c]
    assert len(calls) == 3 * 3                               # 3 groups x 3 chunks
    assert all("--bam-list " in c and "group000" in c and ".bams.txt" in c for c in calls)
    assert not any(" s0.bam" in c for c in calls)            # groups go in via list files
    assert sum(c.startswith("bcftools concat -a ") for c in cmds) == 3   # one per group
    assert sum("bcftools index" in c for c in cmds) == 9
    assert any(c.startswith("# harmonize") for c in cmds)
    assert sum(c.startswith("bcftools annotate") for c in cmds) == 3
    merge = [c for c in cmds if c.startswith("bcftools merge")]
    assert len(merge) == 1 and "AD:sum" in merge[0] and "RPBZ:avg" in merge[0] \
        and "FS:min" in merge[0] and merge[0].count("h_group") == 3
    assert cmds.index(merge[0]) > max(i for i, c in enumerate(cmds) if "annotate" in c)
    assert "reheader" in cmds[-1]                            # naming is still the last step


def test_group_mode_without_regions_is_one_job_per_group(tmp_path):
    cmds = C.call_variants("ref.fa", str(tmp_path / "out.bcf"),
                           bams=[f"s{i}.bam" for i in range(5)], threads=4, bam_batch=2,
                           skip_indels=True, dry_run=True)
    calls = [c for c in cmds if "mpileup" in c]
    assert len(calls) == 3 and not any("-R " in c for c in calls)
    assert not any("concat" in c for c in cmds)              # the group file is the tile
    assert any(c.startswith("bcftools merge") for c in cmds)


def test_a_batch_no_smaller_than_the_cohort_changes_nothing(tmp_path):
    bams = [f"s{i}.bam" for i in range(4)]
    plain = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bams=bams,
                            regions=_bed(tmp_path, 6), threads=2, dry_run=True)
    batched = C.call_variants("ref.fa", str(tmp_path / "o.bcf"), bams=bams,
                              regions=_bed(tmp_path, 6), threads=2, bam_batch=4,
                              dry_run=True)
    # the chunk files live in a fresh temp dir per call; compare everything but that
    strip = lambda cs: [c.split("call_variants.")[0] for c in cs]  # noqa: E731
    assert strip(plain) == strip(batched)
    assert not any("merge" in c for c in batched)


def test_group_mode_needs_the_alignments_listed(tmp_path):
    with pytest.raises(SystemExit, match="--bam-batch needs the alignments listed"):
        C.call_variants("ref.fa", "o.bcf", bam_list=str(tmp_path / "missing.txt"),
                        bam_batch=2, skip_indels=True, dry_run=True)
    with pytest.raises(SystemExit, match="--bam-batch must be 0"):
        C.call_variants("ref.fa", "o.bcf", bams=["a.bam"], bam_batch=-1, dry_run=True)


def test_group_mode_refuses_to_drop_indels_unannounced(tmp_path):
    """Harmonizing is SNP-only. Rather than call indels and quietly lose them, a split
    that is actually in effect demands --skip-indels; a batch the cohort fits in does not,
    since nothing is harmonized and the indels stay."""
    bams = [f"s{i}.bam" for i in range(5)]
    with pytest.raises(SystemExit, match="requires --skip-indels"):
        C.call_variants("ref.fa", "o.bcf", bams=bams, bam_batch=2, dry_run=True)
    cmds = C.call_variants("ref.fa", "o.bcf", bams=bams, bam_batch=2, skip_indels=True,
                           dry_run=True)
    assert all(" -I " in c for c in cmds if "mpileup" in c)
    fits = C.call_variants("ref.fa", "o.bcf", bams=bams, bam_batch=5, dry_run=True)
    assert not any(" -I " in c for c in fits if "mpileup" in c)


def _cohort(tmp_path):
    """Six BAMs whose ALTs are unevenly spread, so a 2-BAM group misses alleles others
    carry: the situation group mode has to get right. Returns (fasta, bams, bed)."""
    pysam = pytest.importorskip("pysam")
    import random
    random.seed(11)
    ref = "".join(random.choice("ACGT") for _ in range(4000))
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
    pysam.faidx(str(fa))
    sites = [1000, 1500, 2000, 2500, 3000]
    alt = lambda b, k: [x for x in "ACGT" if x != b][k]  # noqa: E731
    # sample -> site -> (which alt, fraction of reads carrying it)
    plan = {"s1": {1000: (0, 0.5), 1500: (0, 1.0), 2000: (0, 0.1)},
            "s2": {1000: (0, 0.5), 2500: (1, 1.0)},
            "s3": {1000: (1, 0.5), 3000: (0, 0.5)},
            "s4": {1500: (0, 0.5)},
            "s5": {2500: (0, 0.3), 3000: (0, 1.0)},
            "s6": {}}
    bams = []
    for s, p in plan.items():
        hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
               "SQ": [{"LN": len(ref), "SN": "chr1"}],
               "RG": [{"ID": "rg1", "SM": s, "PL": "ILLUMINA"}]}
        path = tmp_path / f"{s}.bam"
        with pysam.AlignmentFile(str(path), "wb", header=hdr) as out:
            for site in sites:
                for i in range(30):
                    st = site - 80 + i * 2
                    seq = list(ref[st:st + 120])
                    if site in p:
                        k, frac = p[site]
                        if i < int(30 * frac):
                            seq[site - st] = alt(ref[site], k)
                    a = pysam.AlignedSegment()
                    a.query_name = f"r{site}_{i}"; a.query_sequence = "".join(seq)
                    a.flag = 16 if i % 4 in (2, 3) else 0
                    a.reference_id = 0; a.reference_start = st
                    a.mapping_quality = 60 if i % 5 else 40
                    a.cigarstring = "120M"
                    a.query_qualities = pysam.qualitystring_to_array("I" * 120)
                    a.set_tag("RG", "rg1")
                    out.write(a)
        pysam.index(str(path))
        bams.append(str(path))
    bed = tmp_path / "sites.bed"
    bed.write_text("".join(f"chr1\t{s}\t{s + 1}\n" for s in sites))
    return str(fa), bams, str(bed)


def _by_allele(path):
    """Each record keyed by position, with every per-allele value keyed by the allele's
    base rather than its index, so the ALT order (bcftools' by evidence, the harmonized
    union's alphabetical) does not count as a difference."""
    import pysam
    out = {}
    with pysam.VariantFile(path) as vf:
        for r in vf:
            alleles = list(r.alleles)
            per = {"alts": set(r.alts or ()), "DP": r.info["DP"], "DP4": tuple(r.info["DP4"]),
                   "INFO": {f: dict(zip(alleles, r.info[f])) for f in ("AD", "ADF", "ADR")}}
            for s in r.samples:
                gt = r.samples[s]["GT"]
                per[s] = (tuple(sorted(alleles[g] for g in gt)),
                          {f: dict(zip(alleles, r.samples[s][f])) for f in ("AD", "ADF", "ADR")},
                          r.samples[s]["DP"])
            out[(r.chrom, r.pos)] = per
    return out


@needs_bcftools
def test_splitting_the_alignments_gives_the_joint_calls(tmp_path):
    """The property group mode has to have: genotypes, per-sample and per-site allele
    depths and the read counts come out as the joint call had them, for every sample,
    including at sites where a group saw no ALT at all (s6's group at 1001) or a
    different one (C vs G at 2501). Sample order and names are kept; PL is gone."""
    from plasgenomicsutils.lib.bcftools import format_tags, sample_names as names_of
    fa, bams, bed = _cohort(tmp_path)
    joint = str(tmp_path / "joint.bcf")
    grouped = str(tmp_path / "grouped.bcf")
    C.call_variants(fa, joint, bams=bams, regions=bed, threads=2, skip_indels=True)
    C.call_variants(fa, grouped, bams=bams, regions=bed, threads=2, bam_batch=2,
                    skip_indels=True)

    assert names_of(grouped) == names_of(joint) == [f"s{i}" for i in range(1, 7)]
    j, g = _by_allele(joint), _by_allele(grouped)
    assert set(j) == set(g) and len(j) == 5
    assert j == g
    # the union really was needed: no group carried both alleles at 2501, and s6's group
    # had nothing at 1001 -- yet its AD there has a zero for every ALT, not a missing value
    assert g[("chr1", 2501)]["alts"] == {"C", "G"}
    assert g[("chr1", 1001)]["s6"][1]["AD"] == {"A": 30, "C": 0, "G": 0}
    assert "PL" in format_tags(joint) and "PL" not in format_tags(grouped)
    for tag in ("RPBZ", "MQBZ", "FS", "MQ"):          # the by-rule statistics are present
        assert tag in C.info_tags(grouped)


@needs_bcftools
def test_group_mode_is_the_same_with_or_without_region_chunks(tmp_path):
    fa, bams, bed = _cohort(tmp_path)
    a = str(tmp_path / "a.bcf")
    b = str(tmp_path / "b.bcf")
    C.call_variants(fa, a, bams=bams, regions=bed, threads=3, bam_batch=4,
                    skip_indels=True)                                    # 2 x 3 tiles
    C.call_variants(fa, b, bams=bams, regions=bed, threads=1, bam_batch=4,
                    skip_indels=True)                                    # 2 jobs
    assert _by_allele(a) == _by_allele(b)


@needs_bcftools
def test_merge_rules_only_name_tags_the_header_has(tmp_path):
    fa, bams, bed = _cohort(tmp_path)
    out = str(tmp_path / "o.bcf")
    C.call_variants(fa, out, bams=bams[:2], regions=bed)
    rules = C.merge_info_rules(out)
    assert "AD:sum" in rules and "FS:min" in rules and "RPBZ:avg" in rules
    plain = str(tmp_path / "plain.bcf")
    C.call_variants(fa, plain, bams=bams[:2], regions=bed, annotations="FORMAT/AD")
    assert "ADF" not in C.merge_info_rules(plain) and "FS" not in C.merge_info_rules(plain)
    assert "DP:sum" in C.merge_info_rules(plain)
