#!/usr/bin/env python
"""Hard filter on INFO metrics, in GATK's vocabulary or bcftools'; keep PASS."""

from __future__ import annotations

import argparse

from ...lib import vcf_filters as F
from ...lib.bcftools import report_counts
from ...lib.vcf_filters import BCFTOOLS_CALL_RECIPE


def _threshold(v):
    """A threshold, `none` to switch that test off, or `auto` for the caller's default."""
    t = str(v).strip().lower()
    if t == "auto":
        return "auto"
    if t in ("none", "off"):
        return None
    return float(v)


def get_parser_hard_qc_filter() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils hard_qc_filter",
        description="Hard filter on the caller's own QC metrics, keeping PASS. "
                    "--caller gatk (default) reads QD/MQ/SOR/MQRankSum/ReadPosRankSum; "
                    "--caller bcftools reads what bcftools mpileup writes instead "
                    "(FS/RPBZ/SCBZ/MQBZ/MQSBZ), which is the same set of questions asked "
                    "of a bcftools callset.",
        epilog="A bcftools callset needs the right annotations to filter on:\n\n    "
               + BCFTOOLS_CALL_RECIPE + "\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--caller", choices=("gatk", "bcftools"), default="gatk",
                   help="Which caller's metrics the input carries (default: gatk)")
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--qd", type=_threshold, default="auto",
                   help="Drop QD < this (gatk), or QUAL/INFO/DP < this (bcftools). "
                        "'auto' = 20 for gatk and off for bcftools, whose QUAL is not on "
                        "the same scale; 'none' switches it off.")
    p.add_argument("--mq", type=_threshold, default=55.0,
                   help="Drop MQ < this (default: 55; 'none' to switch off)")
    p.add_argument("--sor", type=float, default=3.0, help="Drop SOR > this (default: 3)")
    p.add_argument("--mqranksum", type=float, default=-5.0, help="Drop MQRankSum < this (default: -5)")
    p.add_argument("--readposranksum", type=float, default=-5.0,
                   help="Drop ReadPosRankSum < this (default: -5)")
    p.add_argument("--fs", type=float, default=None,
                   help="gatk only. Optional: drop FS > this (Phred-scaled)")

    b = p.add_argument_group(
        "bcftools metrics (--caller bcftools)",
        "The mirror of the GATK thresholds above. Note two are not straight renames: "
        "bcftools FS is the strand-bias p-value itself, so the test is 'below' a small "
        "number rather than 'above' a large one (Phred 60 is p=1e-6, the same statement); "
        "and the *BZ tags are two-sided, since a read-position artifact leans either way.")
    b.add_argument("--strand-bias-p", type=_threshold, default=1e-6,
                   help="Strand bias: drop INFO/FS below this p-value (default: 1e-6, "
                        "the p-value GATK's FS > 60 corresponds to)")
    b.add_argument("--read-pos-z", type=_threshold, default=5.0,
                   help="Variants sitting at the ends of reads: drop |RPBZ| or |SCBZ| "
                        "above this (default: 5, mirroring ReadPosRankSum < -5)")
    b.add_argument("--max-bias-z", type=_threshold, default=5.0,
                   help="Mapping-quality bias: drop |MQBZ| or |MQSBZ| above this "
                        "(default: 5, mirroring MQRankSum < -5)")
    b.add_argument("--bqbz-z", type=_threshold, default=None,
                   help="Optional: drop |BQBZ| (base-quality bias) above this")
    b.add_argument("--mq0f", type=_threshold, default=None,
                   help="Optional: drop MQ0F (fraction of MQ0 reads) above this")
    p.add_argument("--keep-bed", default=None,
                   help="Whitelist BED of regions to keep whatever this filter says "
                        "(0-based half-open, like any BED). Whitelisted variants still "
                        "face every other filter -- this only exempts them from this "
                        "one.")
    return p


def parse_args_hard_qc_filter():
    return get_parser_hard_qc_filter().parse_args()


def hard_qc_filter():
    args = parse_args_hard_qc_filter()
    F.hard_qc_filter(args.input, args.output, caller=args.caller, qd=args.qd, mq=args.mq,
                     sor=args.sor, mqranksum=args.mqranksum,
                     readposranksum=args.readposranksum, fs=args.fs,
                     strand_bias_p=args.strand_bias_p, read_pos_z=args.read_pos_z,
                     max_bias_z=args.max_bias_z, bqbz_z=args.bqbz_z, mq0f=args.mq0f,
                     keep_bed=args.keep_bed)
    report_counts(args.input, args.output, "hard_qc_filter")


if __name__ == "__main__":
    hard_qc_filter()
