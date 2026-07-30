"""Single runner / dispatch entrypoint.

One CLI, subcommands grouped by area. ``REGISTRY`` maps
group -> command_name -> Command(func, help). Commands share a single flat,
globally-unique namespace on the command line; the group is shown in ``--help``
and in the grouped ``--list`` catalog.

To add a command: write a leaf under ``scripts/<group>/`` exposing a zero-arg
entry function, import it here, and add one ``REGISTRY`` line.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from typing import Callable, Dict

from . import __version__

# -- IBD leaves ---------------------------------------------------------------
from .scripts.ibd.build_matrix import build_matrix
from .scripts.ibd.compute_allele_freqs import compute_allele_freqs
from .scripts.ibd.analyze_matrix import analyze_matrix
from .scripts.ibd.selection_statistic import selection_statistic
from .scripts.ibd.fraction_and_snp_density import fraction_and_snp_density

# -- VCF leaves ---------------------------------------------------------------
from .scripts.vcf.filter_ad_regenotype import filter_ad_regenotype
from .scripts.vcf.harmonize_bcf import harmonize_bcf
from .scripts.vcf.hard_qc_filter import hard_qc_filter
from .scripts.vcf.singleton_filter_add_ads import singleton_filter_add_ads
from .scripts.vcf.biallelic_snp_filter import biallelic_snp_filter
from .scripts.vcf.tandem_repeat_mask import tandem_repeat_mask
from .scripts.vcf.core_region_filter import core_region_filter
from .scripts.vcf.paralog_mask import paralog_mask
from .scripts.vcf.sample_coverage_filter import sample_coverage_filter
from .scripts.vcf.locus_missingness_filter import locus_missingness_filter
from .scripts.vcf.maf_filter import maf_filter
from .scripts.vcf.filter_pipeline import filter_pipeline
from .scripts.vcf.strand_bias_scan import strand_bias_scan
from .scripts.vcf.strand_read_check import strand_read_check


@dataclass(frozen=True)
class Command:
    func: Callable[[], None]  # zero-arg; the leaf parses its own args
    help: str


# group -> command_name -> Command.  Command names must be globally unique.
REGISTRY: Dict[str, Dict[str, Command]] = {
    "ibd": {
        "build_ibd_matrix": Command(build_matrix,
            "Build binary (pairs x SNPs) IBD matrix from hmmibd-rs blocks"),
        "compute_allele_freqs": Command(compute_allele_freqs,
            "Compute global + per-region allele frequencies (single pass)"),
        "analyze_ibd_matrix": Command(analyze_matrix,
            "Per-pair / per-SNP / per-region / per-chr IBD summaries"),
        "ibd_selection_statistic": Command(selection_statistic,
            "IBD-based selection statistic (XiR,s), genome-wide and per-region"),
        "ibd_fraction_and_snp_density": Command(fraction_and_snp_density,
            "Per-pair IBD fraction (callable denominator) and SNP density"),
    },
    "vcf": {
        "filter_ad_regenotype": Command(filter_ad_regenotype,
            "Clean within-sample AD artifacts by depth/frequency, then re-genotype"),
        "harmonize_bcf": Command(harmonize_bcf,
            "Harmonize ALT sets of separately-called cohorts for bcftools merge"),
        "hard_qc_filter": Command(hard_qc_filter,
            "GATK-style hard filter on INFO metrics (QD/MQ/SOR/RankSums), keep PASS"),
        "singleton_filter_add_ads": Command(singleton_filter_add_ads,
            "Drop near-private variants and add the FORMAT/ADS summed-depth tag"),
        "biallelic_snp_filter": Command(biallelic_snp_filter,
            "Keep biallelic SNPs, trimming ALT alleles unused after re-genotyping"),
        "tandem_repeat_mask": Command(tandem_repeat_mask,
            "Remove variants overlapping a tandem-repeat BED"),
        "core_region_filter": Command(core_region_filter,
            "Keep only variants inside the core-genome BED"),
        "paralog_mask": Command(paralog_mask,
            "Remove variants overlapping paralogous/multigene-family genes"),
        "sample_coverage_filter": Command(sample_coverage_filter,
            "Drop low-coverage samples; refresh AC/AN/AF"),
        "locus_missingness_filter": Command(locus_missingness_filter,
            "Keep loci with low missingness and high per-sample coverage"),
        "maf_filter": Command(maf_filter,
            "Keep variants within a minor-allele-frequency window"),
        "filter_pipeline": Command(filter_pipeline,
            "Run an ordered, config-driven chain of filtering steps, tallying counts"),
        "strand_bias_scan": Command(strand_bias_scan,
            "Flag strand-bias (SSE) fake-het artifacts from FORMAT/ADF+ADR; emit a blacklist BED"),
        "strand_read_check": Command(strand_read_check,
            "Read-level strand-bias diagnostic at one site (+ optional ALT-read extraction)"),
    },
}


def _flatten() -> dict[str, tuple[str, Command]]:
    index: dict[str, tuple[str, Command]] = {}
    for group, commands in REGISTRY.items():
        for name, cmd in commands.items():
            if name in index:
                raise RuntimeError(f"Duplicate command name detected: '{name}'")
            index[name] = (group, cmd)
    return index


def _print_catalog() -> None:
    print(f"plasgenomicsutils {__version__} — Plasmodium genomics post-processing\n")
    print("Usage: plasgenomicsutils <command> [options]")
    print("       plasgenomicsutils --list          # this catalog")
    print("       plasgenomicsutils <command> -h    # help for one command\n")
    for group, commands in REGISTRY.items():
        print(f"[{group}]")
        width = max((len(n) for n in commands), default=0)
        for name, cmd in commands.items():
            print(f"  {name.ljust(width)}   {cmd.help}")
        print()


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    index = _flatten()

    if not argv or argv[0] in ("-h", "--help", "--list"):
        _print_catalog()
        return
    if argv[0] in ("-V", "--version"):
        print(f"plasgenomicsutils {__version__}")
        return

    command = argv[0]
    if command not in index:
        print(f"ERROR: unknown command '{command}'\n", file=sys.stderr)
        _print_catalog()
        raise SystemExit(2)

    _group, cmd = index[command]
    leaf_prog = f"plasgenomicsutils {command}"
    old_argv = sys.argv[:]
    try:
        sys.argv = [leaf_prog, *argv[1:]]  # leaf parses its own args
        cmd.func()
    finally:
        sys.argv = old_argv


if __name__ == "__main__":
    main()
