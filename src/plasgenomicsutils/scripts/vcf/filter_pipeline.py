#!/usr/bin/env python
"""Run an ordered, config-driven chain of VCF filtering steps, tallying counts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ...lib import filter_pipeline as P


def get_parser_filter_pipeline() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils filter_pipeline",
        description="Run an ordered chain of filtering steps defined by a JSON config.",
    )
    p.add_argument("--input", help="Input VCF/BCF")
    p.add_argument("--config", help="Pipeline JSON config (see --emit-default-config)")
    p.add_argument("--outdir", default="filter_pipeline", help="Output directory")
    p.add_argument("--emit-default-config", metavar="PATH", default=None,
                   help="Write a starting default config to PATH and exit")
    p.add_argument("--no-snp-bed", action="store_true",
                   help="Do not auto-write the final SNP-panel BED (used by the IBD tools)")
    return p


def parse_args_filter_pipeline():
    return get_parser_filter_pipeline().parse_args()


def _tally_fields(row: dict) -> tuple[str, "int | str", str]:
    """``(kind, count, path)`` for one ``run_pipeline()`` tally row.

    The rows are not all one shape: a filter step records ``variants``, a report records
    ``rows``, and a step switched off with ``"enabled": false`` records neither.
    """
    if row.get("skipped"):
        return "skipped", "", ""
    if row.get("report"):
        return "report", row["rows"], row.get("path", "")
    return "variants", row["variants"], row.get("path", "")


def filter_pipeline():
    args = parse_args_filter_pipeline()

    if args.emit_default_config:
        Path(args.emit_default_config).write_text(json.dumps(P.DEFAULT_CONFIG, indent=2) + "\n")
        print(f"Wrote default config -> {args.emit_default_config}")
        return

    if not args.input or not args.config:
        raise SystemExit("ERROR: --input and --config are required (or use --emit-default-config)")

    config = P.load_config(args.config)
    tally = P.run_pipeline(args.input, args.outdir, config, emit_snp_bed=not args.no_snp_bed)

    print("\n=== variant counts per step ===")
    for row in tally:
        kind, count, _ = _tally_fields(row)
        shown = "skipped" if kind == "skipped" else f"{count:,}"
        print(f"  {row['step']:<28} {shown:>12}{' rows' if kind == 'report' else ''}")

    lines = ["step\tkind\tcount\tpath\n"]
    for row in tally:
        kind, count, path = _tally_fields(row)
        lines.append(f"{row['step']}\t{kind}\t{count}\t{path}\n")
    Path(args.outdir, "variant_counts.tsv").write_text("".join(lines))
    print(f"\n  -> {Path(args.outdir, 'variant_counts.tsv')}")


if __name__ == "__main__":
    filter_pipeline()
