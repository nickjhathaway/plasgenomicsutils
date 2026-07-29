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
    # "vcf": { ... }
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
