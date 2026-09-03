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
from .scripts.ibd.gene_overlap import gene_overlap
from .scripts.ibd.gene_pairs import gene_pairs

# -- Fws leaves ---------------------------------------------------------------
from .scripts.fws.calculate_fws import calculate_fws

# -- LD leaves ----------------------------------------------------------------
from .scripts.ld.decay import ld_decay

# -- Coverage leaves ----------------------------------------------------------
from .scripts.cov.depth_stats import depth_stats
from .scripts.cov.dropout_regions import dropout_regions

# -- VCF leaves ---------------------------------------------------------------
from .scripts.vcf.call_variants import call_variants
from .scripts.vcf.filter_ad_regenotype import filter_ad_regenotype
from .scripts.vcf.fws_filter import fws_filter
from .scripts.vcf.harmonize_bcf import harmonize_bcf
from .scripts.vcf.no_alt_filter import no_alt_filter
from .scripts.vcf.hard_qc_filter import hard_qc_filter
from .scripts.vcf.singleton_filter_add_ads import singleton_filter_add_ads
from .scripts.vcf.singleton_counts import singleton_counts
from .scripts.vcf.wsaf_profile import wsaf_profile
from .scripts.vcf.biallelic_snp_filter import biallelic_snp_filter
from .scripts.vcf.strip_stale_format import strip_stale_format
from .scripts.vcf.tandem_repeat_mask import tandem_repeat_mask
from .scripts.vcf.core_region_filter import core_region_filter
from .scripts.vcf.paralog_mask import paralog_mask
from .scripts.vcf.sample_coverage_filter import sample_coverage_filter
from .scripts.vcf.locus_missingness_filter import locus_missingness_filter
from .scripts.vcf.maf_filter import maf_filter
from .scripts.vcf.filter_pipeline import filter_pipeline
from .scripts.vcf.strand_bias_scan import strand_bias_scan
from .scripts.vcf.variant_spacing import variant_spacing
from .scripts.vcf.vcf_to_bed import vcf_to_bed
from .scripts.vcf.strand_read_check import strand_read_check


@dataclass(frozen=True)
class Command:
    func: Callable[[], None]  # zero-arg; the leaf parses its own args
    help: str


# group -> command_name -> Command.  Command names must be globally unique.
REGISTRY: Dict[str, Dict[str, Command]] = {
    "fws": {
        "calculate_fws": Command(calculate_fws,
            "Per-sample Fws within-host diversity (moimix::getFws) from a VCF or AD table"),
    },
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
        "ibd_gene_overlap": Command(gene_overlap,
            "Fraction of pairs whose IBD block overlaps each gene, per group pair"),
        "ibd_gene_pairs": Command(gene_pairs,
            "Sample pairs IBD over each gene, with how much of the gene is covered"),
    },
    "ld": {
        "ld_decay": Command(ld_decay,
            "Mean r-squared vs SNP-pair distance per group: how fast LD decays"),
    },
    "coverage": {
        "coverage_depth_stats": Command(depth_stats,
            "Per-sample depth: mean/median/SD and breadth at thresholds, per chromosome"),
        "coverage_dropout_regions": Command(dropout_regions,
            "Regions below depth in nearly every sample (sWGA amplification dropouts)"),
    },
    # Working with a callset itself: making one, merging one, converting one.
    "vcf": {
        "call_variants": Command(call_variants,
            "Call variants with bcftools, annotated for hard_qc_filter, parallel over regions"),
        "harmonize_bcf": Command(harmonize_bcf,
            "Harmonize ALT sets of separately-called cohorts for bcftools merge"),
        "vcf_to_bed": Command(vcf_to_bed,
            "Convert a VCF/BCF to 0-based BED (stdout by default)"),
    },
    # The filtering chain: the runner, then its steps **in the order the default
    # config runs them** rather than alphabetically -- the order is the point, and a
    # step reads differently depending on what has already happened. Every one is also
    # a standalone command.
    "vcf_filter_pipeline": {
        "filter_pipeline": Command(filter_pipeline,
            "Run an ordered, config-driven chain of filtering steps, tallying counts"),
        "no_alt_filter": Command(no_alt_filter,
            "Drop records with no ALT allele (non-variant positions), counted separately"),
        "hard_qc_filter": Command(hard_qc_filter,
            "GATK-style hard filter on INFO metrics (QD/MQ/SOR/RankSums), keep PASS"),
        "singleton_filter_add_ads": Command(singleton_filter_add_ads,
            "Drop near-private variants and add the FORMAT/ADS summed-depth tag"),
        "tandem_repeat_mask": Command(tandem_repeat_mask,
            "Remove variants overlapping a tandem-repeat BED"),
        "core_region_filter": Command(core_region_filter,
            "Keep only variants inside the core-genome BED"),
        "paralog_mask": Command(paralog_mask,
            "Remove variants overlapping paralogous/multigene-family genes"),
        "filter_ad_regenotype": Command(filter_ad_regenotype,
            "Clean within-sample AD artifacts by depth/frequency, then re-genotype"),
        "biallelic_snp_filter": Command(biallelic_snp_filter,
            "Keep biallelic SNPs, trimming ALT alleles unused after re-genotyping"),
        "sample_coverage_filter": Command(sample_coverage_filter,
            "Drop low-coverage samples; refresh AC/AN/AF"),
        "locus_missingness_filter": Command(locus_missingness_filter,
            "Keep loci with low missingness and high per-sample coverage"),
        "maf_filter": Command(maf_filter,
            "Keep variants within a minor-allele-frequency window"),
        "fws_filter": Command(fws_filter,
            "Keep only monoclonal samples (Fws >= a threshold); refresh AC/AN/AF"),
        "strip_stale_format": Command(strip_stale_format,
            "Strip stale genotype-linked FORMAT fields (e.g. PL) that no longer match the genotypes"),
    },
    # Reads a callset and writes a table; never changes the data.
    "vcf_reporting": {
        "singleton_counts": Command(singleton_counts,
            "Per-sample count of variants where it is the only non-reference carrier"),
        "strand_bias_scan": Command(strand_bias_scan,
            "Flag strand-bias (SSE) fake-het artifacts from FORMAT/ADF+ADR; emit a blacklist BED"),
        "strand_read_check": Command(strand_read_check,
            "Read-level strand-bias diagnostic at one site (+ optional ALT-read extraction)"),
        "variant_spacing": Command(variant_spacing,
            "Per-chromosome gaps between consecutive variants, and density per cM"),
        "wsaf_profile": Command(wsaf_profile,
            "Per-sample: is there a dominant clone, and what --min-freq reduces the sample to it"),
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


def _print_catalog_plain() -> None:
    """Machine-friendly '<command>\t<group>\t<help>' (one per line), for bash completion."""
    for group, commands in REGISTRY.items():
        for name, cmd in commands.items():
            print(f"{name}\t{group}\t{cmd.help}")


# Canonical bash-completion script. `plasgenomicsutils --bash-completion` prints this, and
# etc/bash_completion in the repo is a copy for people who'd rather source a file directly.
_BASH_COMPLETION = r"""# bash completion for plasgenomicsutils
# enable with:
#   plasgenomicsutils --bash-completion >> ~/.bash_completion && source ~/.bash_completion

_plasgenomicsutils_complete()
{
    # keep '_' and '=' out of the word breaks so --group-col / --bcf=path stay one token
    local _OLD_WB="${COMP_WORDBREAKS-}"
    COMP_WORDBREAKS="${COMP_WORDBREAKS//_/}"
    COMP_WORDBREAKS="${COMP_WORDBREAKS//=}"

    local cur
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"

    # 1) first token: the command name (queried live from the CLI)
    if [[ ${COMP_CWORD} -eq 1 ]]; then
        local cmds
        cmds="$("${COMP_WORDS[0]}" --list-plain 2>/dev/null | awk -F'\t' '{print $1}')"
        COMPREPLY=( $(compgen -W "${cmds}" -- "${cur}") )
        COMP_WORDBREAKS="$_OLD_WB"
        return 0
    fi

    # 2) an option for a leaf command: scrape its -h output for -x / --long flags.
    # Keep only the flag portion of each option line (strip the indent, then cut at the
    # 2+ spaces before the help text) so words like "hmmibd-rs" in a description can't
    # masquerade as flags.
    if [[ "${cur}" == -* ]]; then
        local opts
        opts="$("${COMP_WORDS[0]}" "${COMP_WORDS[1]}" -h 2>/dev/null \
            | awk '/^[[:space:]]+-/ { sub(/^[[:space:]]+/, ""); sub(/  .*/, ""); print }' \
            | grep -oE -- '-{1,2}[A-Za-z0-9][-A-Za-z0-9_]*' | sort -u)"
        COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
        COMP_WORDBREAKS="$_OLD_WB"
        return 0
    fi

    # 3) otherwise, filename completion for positional args
    COMPREPLY=( $(compgen -f -- "${cur}") )
    COMP_WORDBREAKS="$_OLD_WB"
    return 0
}

complete -F _plasgenomicsutils_complete plasgenomicsutils
"""


def _print_bash_completion() -> None:
    sys.stdout.write(_BASH_COMPLETION)


def main(argv: list[str] | None = None) -> None:
    """Dispatch one subcommand, or print the grouped catalog for ``--list``."""
    argv = list(sys.argv[1:] if argv is None else argv)
    index = _flatten()

    if not argv or argv[0] in ("-h", "--help", "--list"):
        _print_catalog()
        return
    if argv[0] == "--list-plain":
        _print_catalog_plain()
        return
    if argv[0] == "--bash-completion":
        _print_bash_completion()
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
