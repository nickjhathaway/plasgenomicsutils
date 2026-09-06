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
                    "(RPBZ/SCBZ/MQBZ/MQSBZ, and SOR computed from ADF/ADR), which is the "
                    "same set of questions asked of a bcftools callset.",
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
    p.add_argument("--sor", type=_threshold, default=3.0,
                   help="Strand bias, in both modes: drop SOR > this (default: 3, GATK's "
                        "own cutoff). GATK writes SOR; for a bcftools callset it is "
                        "computed from INFO/ADF and INFO/ADR, so the threshold means the "
                        "same thing either way. 'none' switches it off.")
    p.add_argument("--mqranksum", type=float, default=-5.0, help="Drop MQRankSum < this (default: -5)")
    p.add_argument("--readposranksum", type=float, default=-5.0,
                   help="Drop ReadPosRankSum < this (default: -5)")
    p.add_argument("--fs", type=float, default=None,
                   help="gatk only. Optional: drop FS > this (Phred-scaled)")

    b = p.add_argument_group(
        "bcftools metrics (--caller bcftools)",
        "The mirror of the GATK thresholds above. Strand bias is --sor in both modes; "
        "these are the rest. Note the *BZ tags are two-sided, since a read-position "
        "artifact leans either way, and that they are test statistics: a z-score grows "
        "with the reads pooled across samples, so a cutoff that suits 20 samples is "
        "stricter at 400. Each z test therefore also asks how big the shift behind the z "
        "is (--bias-eff), which does not move with depth.")
    b.add_argument("--strand-bias-p", type=_threshold, default="auto",
                   help="Strand bias as a significance test instead: drop INFO/FS below "
                        "this p-value. Off by default ('auto'), because FS is a p-value "
                        "pooled over every sample -- at a few hundred samples it rejects "
                        "sites with no meaningful skew, and under ~1e-38 it saturates at "
                        "0.0 in the file's 32-bit float. Use --sor. 1e-6 is what GATK's "
                        "FS > 60 corresponds to, if you want it anyway.")
    b.add_argument("--read-pos-z", type=_threshold, default=5.0,
                   help="Variants sitting at the ends of reads: drop |RPBZ| or |SCBZ| "
                        "above this (default: 5, mirroring ReadPosRankSum < -5), provided "
                        "the shift is also at least --bias-eff")
    b.add_argument("--max-bias-z", type=_threshold, default=5.0,
                   help="Mapping-quality bias: drop |MQBZ| or |MQSBZ| above this "
                        "(default: 5, mirroring MQRankSum < -5)")
    b.add_argument("--bqbz-z", type=_threshold, default=None,
                   help="Optional: drop |BQBZ| (base-quality bias) above this")
    b.add_argument("--mq0f", type=_threshold, default=None,
                   help="Optional: drop MQ0F (fraction of MQ0 reads) above this")
    b.add_argument("--bias-eff", type=_threshold, default=0.15,
                   help="How big the shift behind a *BZ z-score has to be before the z "
                        "counts, for every z test above. Computed from INFO/ADF+ADR as "
                        "|z|*sqrt((n1+n2+1)/(12*n1*n2)) -- ref vs alt reads for RPBZ/SCBZ/"
                        "MQBZ/BQBZ, forward vs reverse for MQSBZ -- it is how far P(a read "
                        "of one group ranks above a read of the other) sits from 0.5. "
                        "Default 0.15, a 65:35 split (rank-biserial 0.3). This is what "
                        "stops a z test tightening with cohort size: the same shift scores "
                        "z=1 in one sample and z=30 pooled over 400. 'none' restores the "
                        "plain z tests.")
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
                     max_bias_z=args.max_bias_z, bias_eff=args.bias_eff, bqbz_z=args.bqbz_z, mq0f=args.mq0f,
                     keep_bed=args.keep_bed)
    report_counts(args.input, args.output, "hard_qc_filter")


if __name__ == "__main__":
    hard_qc_filter()
