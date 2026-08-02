"""Thin, checked wrappers around bcftools / bedtools for the filtering commands.

External tools (``bcftools``, ``bedtools``, ``bgzip``) must be on ``PATH``. Each
helper validates that the tools it needs are present, runs the command, and
raises a clear error on failure.
"""

from __future__ import annotations

import re
import shlex
import shutil
import subprocess
import sys


def format_tags(path: str) -> set[str]:
    """The FORMAT tag IDs defined in a VCF/BCF header (e.g. ``{"GT", "AD", "PL"}``)."""
    require("bcftools")
    proc = subprocess.run(["bcftools", "view", "-h", str(path)],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return set(re.findall(r"^##FORMAT=<ID=([^,>]+)", proc.stdout, flags=re.M))


def require(*tools: str) -> None:
    """Raise unless every named external tool is on PATH."""
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise SystemExit(
            f"ERROR: required tool(s) not found on PATH: {', '.join(missing)}. "
            "Install bcftools/bedtools (e.g. `conda install -c bioconda bcftools bedtools`)."
        )


def q(path) -> str:
    """Shell-quote a path/argument."""
    return shlex.quote(str(path))


def sh(cmd: str, *, tools: tuple[str, ...] = ()) -> None:
    """Run a shell pipeline (bash), checking required tools first."""
    if tools:
        require(*tools)
    proc = subprocess.run(cmd, shell=True, executable="/bin/bash",
                          stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"ERROR: command failed (exit {proc.returncode}): {cmd}")


def out_flag(path: str) -> str:
    """bcftools -O flag inferred from the output extension."""
    p = str(path)
    if p.endswith(".bcf"):
        return "b"
    if p.endswith(".vcf.gz"):
        return "z"
    return "v"


def index_vcf(path: str) -> None:
    """Create a CSI index for a coordinate-sorted ``.bcf`` / ``.vcf.gz`` (best-effort).

    Indexing intermediate outputs keeps tools that probe for an index (e.g. pysam) quiet
    and lets downstream steps do region queries. Non-fatal: a failure (e.g. an unsorted
    file) is ignored so it never breaks a pipeline.
    """
    p = str(path)
    if not p.endswith((".bcf", ".vcf.gz")):
        return
    try:
        sh(f"bcftools index -f {q(path)}", tools=("bcftools",))
    except SystemExit:
        pass


def count_variants(path: str) -> int:
    """Number of records via ``bcftools view -H | wc -l``."""
    require("bcftools")
    proc = subprocess.run(
        f"bcftools view -H {q(path)} | wc -l",
        shell=True, executable="/bin/bash",
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"ERROR: could not count variants in {path}")
    return int(proc.stdout.strip() or 0)


def report_counts(before: str, after: str, label: str) -> None:
    """Print a per-step ``before -> after`` variant-count line."""
    nb, na = count_variants(before), count_variants(after)
    print(f"  [{label}] variants: {nb:,} -> {na:,}  ({na - nb:+,})")
