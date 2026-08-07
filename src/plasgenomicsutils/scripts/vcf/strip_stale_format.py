#!/usr/bin/env python
"""Strip stale genotype-linked FORMAT fields (default PL) whose stored length no longer
matches the genotypes (e.g. a diploid GT forced over hexaploid calls). Such fields make
`bcftools view --trim-alt-alleles` abort and are meaningless once re-genotyped."""

from __future__ import annotations

import argparse
import sys

from ...lib import strip_format as SF


def get_parser_strip_stale_format() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils strip_stale_format",
        description="Strip stale genotype-linked FORMAT fields (default PL) whose value "
                    "count no longer matches the genotypes, so allele-trimming and "
                    "downstream tools stop choking on them.",
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--fields", nargs="+", default=["PL"], metavar="TAG",
                   help="FORMAT fields to check/strip (default: PL). Genotype-linked "
                        f"candidates: {', '.join(SF.GENOTYPE_LINKED_FORMAT)}.")
    p.add_argument("--mode", choices=["mismatch", "always"], default="mismatch",
                   help="'mismatch' (default): null a field only on records where its "
                        "length is inconsistent with the genotypes (valid records keep "
                        "their values). 'always': drop the fields from every record.")
    return p


def parse_args_strip_stale_format():
    return get_parser_strip_stale_format().parse_args()


def strip_stale_format():
    args = parse_args_strip_stale_format()
    n = SF.strip_stale_format(args.input, args.output,
                              fields=tuple(args.fields), mode=args.mode)
    tags = ", ".join(args.fields)
    if args.mode == "always":
        print(f"[strip_stale_format] dropped {tags} from all records", file=sys.stderr)
    else:
        print(f"[strip_stale_format] nulled stale {tags} on {n} record(s)", file=sys.stderr)


if __name__ == "__main__":
    strip_stale_format()
