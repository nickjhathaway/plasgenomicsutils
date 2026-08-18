#!/usr/bin/env python
"""Per-sample within-sample allele fractions: does a sample have one dominant clone?"""

from __future__ import annotations

import argparse

import pandas as pd

from ...lib.wsaf import (WSAF_MAX_RESIDUAL_HET, WSAF_MIN_DOMINANT, WSAF_MIN_MINOR,
                         WSAF_MIN_MINOR_READS, WSAF_MIN_SITES,
                         wsaf_profile as _profile)
from ...utils.small_utils import Utils


def get_parser_wsaf_profile() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils wsaf_profile",
        description="For each sample, say whether one clone dominates the infection and what "
                    "filter_ad_regenotype --min-freq would reduce it to that clone. Fws says "
                    "how clonal a sample is; this says whether a sample failing that gate can "
                    "still be used. Not a COI estimate.",
    )
    p.add_argument("--input-vcf", required=True, help="VCF/BCF with FORMAT/AD")
    p.add_argument("--out", default="-",
                   help="Per-sample summary TSV ('-' = STDOUT, default)")
    p.add_argument("--sites-out", default=None,
                   help="Also write the per-(sample, site) fractions behind the summary, "
                        "for plotting the distributions (plot_wsaf() in plasgenomicsutilsR)")
    p.add_argument("--fws", default=None,
                   help="Optional calculate_fws output. Focuses the report on the samples "
                        "that fail the gate, which is the question worth asking: of those, "
                        "how many have a dominant clone?")
    p.add_argument("--monoclonal-threshold", type=float, default=0.92,
                   help="Fws at or above which a sample already counts as monoclonal, used "
                        "with --fws (default: 0.92)")
    p.add_argument("--min-dominant", type=float, default=WSAF_MIN_DOMINANT,
                   help=f"Share of the parasitaemia the dominant clone must hold to be worth "
                        f"re-genotyping to (default: {WSAF_MIN_DOMINANT}). The complement is "
                        f"the filter: 0.70 here means --min-freq 0.30.")
    p.add_argument("--max-residual-het", type=float, default=WSAF_MAX_RESIDUAL_HET,
                   help=f"Fraction of covered sites allowed to stay heterozygous after that "
                        f"filter (default: {WSAF_MAX_RESIDUAL_HET})")
    p.add_argument("--min-depth", type=int, default=10,
                   help="Read depth a site needs before its fractions are used (default: "
                        "10; below this the fraction is too coarsely quantised to place a "
                        "band)")
    p.add_argument("--min-minor", type=float, default=WSAF_MIN_MINOR,
                   help=f"Minor fraction a site needs to count as heterozygous (default: "
                        f"{WSAF_MIN_MINOR}, which excludes the stray reads ordinary error "
                        f"rates produce)")
    p.add_argument("--min-minor-reads", type=int, default=WSAF_MIN_MINOR_READS,
                   help=f"Reads the minor allele needs as well as --min-minor (default: "
                        f"{WSAF_MIN_MINOR_READS}). Matters at low depth, where one stray "
                        f"read is 10%% of a 10x site.")
    p.add_argument("--min-sites", type=int, default=WSAF_MIN_SITES,
                   help=f"Covered sites needed to judge a sample rather than call it "
                        f"undetermined (default: {WSAF_MIN_SITES})")
    return p


def parse_args_wsaf_profile():
    return get_parser_wsaf_profile().parse_args()


def _report(df: pd.DataFrame, min_dominant: float) -> None:
    counts = df["class"].value_counts()
    for k in ("monoclonal", "dominant_clone", "mixed", "undetermined"):
        if k in counts:
            print(f"  {k:<16} {counts[k]:>5}")
    dom = df[df["class"] == "dominant_clone"]
    if len(dom):
        need = dom["min_freq_needed"].max()
        print(f"\n  {len(dom)} sample(s) have a dominant clone holding {min_dominant:.0%} or "
              f"more. One filter at\n  --min-freq {need:.2f} covers all of them; per sample "
              f"the threshold it needs is its own\n  min_freq_needed. Re-run calculate_fws on "
              f"the filtered callset and gate on that.")
    mixed = df[df["class"] == "mixed"]
    if len(mixed):
        print(f"\n  {len(mixed)} sample(s) carry a strain too large to remove: no --min-freq "
              f"below 0.5\n  reduces them to one clone without deleting a real strain.")


def wsaf_profile():
    args = parse_args_wsaf_profile()
    df = _profile(args.input_vcf, min_depth=args.min_depth, min_minor=args.min_minor,
                  min_minor_reads=args.min_minor_reads, min_dominant=args.min_dominant,
                  max_residual_het=args.max_residual_het, min_sites=args.min_sites,
                  sites_out=args.sites_out)

    if args.fws:
        fws = Utils.read_table(args.fws)
        if not {"sample", "fws"} <= set(fws.columns):
            raise SystemExit(f"{args.fws} needs `sample` and `fws` columns; found: "
                             f"{', '.join(map(str, fws.columns))}")
        fws["sample"] = fws["sample"].astype(str)
        df["sample"] = df["sample"].astype(str)
        df = df.merge(fws[["sample", "fws"]], on="sample", how="left")
        failing = df[df["fws"] < args.monoclonal_threshold]
        print(f"  {len(failing)} of {len(df)} samples fail Fws >= "
              f"{args.monoclonal_threshold:.2f}. Of those:\n")
        _report(failing, args.min_dominant)
        print(f"\n  (the other {len(df) - len(failing)} already pass the gate; the table "
              f"written out covers every sample)")
    else:
        _report(df, args.min_dominant)

    with Utils.smart_open_write(args.out) as fh:
        df.to_csv(fh, sep="\t", index=False)
    if args.out != "-":
        print(f"\n  -> {args.out}")
    if args.sites_out:
        print(f"  -> {args.sites_out}")


if __name__ == "__main__":
    wsaf_profile()
