#!/usr/bin/env python
"""Call variants with bcftools, annotated for hard_qc_filter, parallel over regions."""

from __future__ import annotations

import argparse

from ...lib import call_variants as C
from ...lib.bcftools import count_variants, sample_names
from ...lib.vcf_filters import BCFTOOLS_MPILEUP_ANNOTATIONS


def get_parser_call_variants() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils call_variants",
        description="Run bcftools mpileup | bcftools call, asking for the annotations "
                    "hard_qc_filter --caller bcftools reads, and splitting a region list "
                    "over concurrent jobs. Calling is one core per job, so parallelism "
                    "comes from splitting the regions rather than from --threads inside "
                    "bcftools.",
        epilog="Example -- 400 SNPs of interest over 8 cores:\n\n"
               "    plasgenomicsutils call_variants --ref Pf3D7.fasta \\\n"
               "      --bam-list bams.txt --regions crt_region_snps.bed \\\n"
               "      --threads 8 --output crt_snps.bcf\n\n"
               "The region list is split into 8 chunks, one job each; the parts are "
               "concatenated and indexed, and\nthe samples are named after their BAM "
               "files rather than their paths.\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--ref", required=True, help="Reference FASTA (indexed)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--bam", nargs="+", help="One or more BAM/CRAM files")
    src.add_argument("--bam-list", help="File of BAM paths, one per line")
    p.add_argument("--output", required=True, help="Output .bcf / .vcf.gz")

    r = p.add_argument_group("regions and parallelism")
    r.add_argument("--regions", default=None,
                   help="Region file to call over. A '.bed' is read 0-based half-open "
                        "and anything else as 1-based CHROM POS -- bcftools' own rule, "
                        "kept when the file is split.")
    r.add_argument("--threads", type=int, default=1,
                   help="Concurrent jobs over chunks of --regions (default: 1). This is "
                        "what makes calling faster; bcftools' own --threads only "
                        "parallelises compression.")
    r.add_argument("--chunk-size", type=int, default=None,
                   help="Regions per chunk. Default splits the list into --threads "
                        "pieces; set this for fixed-size chunks instead, which evens out "
                        "uneven regions at the cost of more jobs.")
    r.add_argument("--keep-chunks", default=None,
                   help="Directory to keep the per-chunk region files and BCFs in "
                        "(default: a temporary directory, removed after concatenating)")

    c = p.add_argument_group("calling")
    c.add_argument("--ploidy", default="2",
                   help="bcftools call --ploidy (default: 2, matching the conventional "
                        "diploid coding the rest of the package expects)")
    c.add_argument("--ignore-rg", action=argparse.BooleanOptionalAction, default=True,
                   help="Treat each alignment as one sample whatever its read groups say "
                        "(default: on, which is the usual one-BAM-per-sample case). "
                        "--no-ignore-rg reads sample names from the RG SM tags instead.")
    c.add_argument("--sample-suffix", default=".bam",
                   help="With --ignore-rg, bcftools names each sample after its path; the "
                        "samples are renamed to the file name with this removed "
                        "(default: '.bam'). Pass the whole trailing part you want gone, "
                        "e.g. --sample-suffix .sorted.dup.pf.bam")
    c.add_argument("--no-rename-samples", action="store_true",
                   help="Keep the path-derived sample names bcftools writes, instead of "
                        "renaming them")
    c.add_argument("--variants-only", action="store_true",
                   help="bcftools call -v: emit only variant sites. Off by default, since "
                        "calling a list of known positions is usually about filling them "
                        "in -- a reference call at a target position is the answer 'this "
                        "sample is reference here', and -v would drop it.")
    c.add_argument("--skip-indels", action="store_true",
                   help="mpileup -I: do not call indels at all, so no INDEL records are "
                        "produced. SNP records are unaffected -- reads spanning an indel "
                        "still count toward the pileup at every other position -- so this "
                        "is not the same as filtering indels out afterwards only in that "
                        "the indel likelihoods are never computed.")
    c.add_argument("--max-depth", type=int, default=None,
                   help="mpileup -d: per-file depth cap for the pileup (bcftools default "
                        "250). Raise it for deep data, or positions are downsampled.")
    c.add_argument("--min-mapq", type=int, default=None,
                   help="mpileup -q: skip alignments with mapping quality below this")
    c.add_argument("--min-baseq", type=int, default=None,
                   help="mpileup -Q: skip bases with base quality below this")
    c.add_argument("--annotations", default=BCFTOOLS_MPILEUP_ANNOTATIONS,
                   help="mpileup -a list (default: what hard_qc_filter --caller bcftools "
                        "reads)")
    c.add_argument("--extra-mpileup", default="",
                   help="Extra arguments passed through to bcftools mpileup")
    c.add_argument("--extra-call", default="",
                   help="Extra arguments passed through to bcftools call")

    p.add_argument("--dry-run", action="store_true",
                   help="Print the commands that would run, and stop")
    return p


def parse_args_call_variants():
    return get_parser_call_variants().parse_args()


def call_variants():
    args = parse_args_call_variants()
    if args.regions and args.threads > 1:
        n = C.n_regions(args.regions)
        size = args.chunk_size or -(-n // args.threads)
        print(f"  {n} region(s) over {args.threads} job(s), ~{size} per chunk")

    cmds = C.call_variants(
        args.ref, args.output, bams=args.bam, bam_list=args.bam_list,
        regions=args.regions, threads=args.threads, chunk_size=args.chunk_size,
        annotations=args.annotations, ploidy=args.ploidy, ignore_rg=args.ignore_rg,
        sample_suffix=None if args.no_rename_samples else args.sample_suffix,
        skip_indels=args.skip_indels, max_depth=args.max_depth, min_mapq=args.min_mapq,
        min_baseq=args.min_baseq, variants_only=args.variants_only,
        extra_mpileup=args.extra_mpileup,
        extra_call=args.extra_call, keep_chunks=args.keep_chunks, dry_run=args.dry_run)

    if args.dry_run:
        for c in cmds:
            print(c)
        return
    names = sample_names(args.output)
    shown = ", ".join(names[:4]) + (f", ... (+{len(names) - 4})" if len(names) > 4 else "")
    print(f"  {count_variants(args.output)} variants, {len(names)} sample(s): {shown}")
    print(f"  -> {args.output}")


if __name__ == "__main__":
    call_variants()
