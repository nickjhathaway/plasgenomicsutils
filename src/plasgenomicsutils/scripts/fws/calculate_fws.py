#!/usr/bin/env python
"""Compute the Fws within-host diversity statistic from a VCF/BCF or an AD table."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ...lib.fws import (MULTIALLELIC_MODES, compute_fws, load_exclude_regions,
                        read_ad_table, read_ad_vcf, read_allele_table, read_pop_freqs,
                        write_pop_freqs)


def get_parser_calculate_fws() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils calculate_fws",
        description="Per-sample Fws (Manske 2012 / moimix::getFws) from per-sample allele "
                    "depths, read from a VCF/BCF, a bcftools-query AD table, or a "
                    "long-format amplicon allele table (microhaplotypes).")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--input-vcf", help="VCF/BCF with FORMAT/AD. Pass it unsplit: a "
                                         "multiallelic record is collapsed with every "
                                         "allele's depth counted (see --multiallelic), "
                                         "whereas `bcftools norm -m-` output has already "
                                         "thrown away the other alleles' reads")
    src.add_argument("--ad-table", help="bcftools query TSV: CHROM POS REF ALT then one full "
                                        "AD ('ref,alt[,alt2..]') per sample")
    src.add_argument("--allele-table", help="Long-format TSV with one row per sample, locus "
                                            "and allele (microhaplotypes): each locus is "
                                            "scored as one multiallelic site whose alleles "
                                            "are its haplotypes. Column names default to "
                                            "MAD4HATTER's; see --sample-col etc. Pair with "
                                            "--n-bins 0.")
    p.add_argument("--samples", help="For --ad-table: comma-separated sample ids in column "
                                     "order, or a file with one per line")
    p.add_argument("--sample-col", default="library_sample_name",
                   help="--allele-table: sample column (default: library_sample_name)")
    p.add_argument("--locus-col", default="target_name",
                   help="--allele-table: locus column (default: target_name)")
    p.add_argument("--allele-col", default="seq",
                   help="--allele-table: allele/haplotype column (default: seq)")
    p.add_argument("--reads-col", default="reads",
                   help="--allele-table: read-count column (default: reads)")
    p.add_argument("--estimator", choices=["regression", "ratio"], default="regression",
                   help="regression = moimix::getFws (default); ratio = summed binned-mean ratio")
    p.add_argument("--min-depth", type=int, default=0,
                   help="Drop per-sample sites below this read depth (default: 0; CNV gate uses 10)")
    p.add_argument("--pop-freqs", help="Use these population allele frequencies instead of "
                                       "the input's own pooled read fractions: a TSV with "
                                       "one row per locus and allele (columns per "
                                       "--freq-*-col). Locus is CHROM:POS for VCF/AD-table "
                                       "input, the locus name for --allele-table; the allele "
                                       "is the REF/ALT string or the haplotype. Sites with no "
                                       "entry are dropped. --write-pop-freqs makes one.")
    p.add_argument("--write-pop-freqs", help="Write the input's own population allele "
                                             "frequencies to this TSV, for reuse as "
                                             "--pop-freqs on another callset")
    p.add_argument("--freq-locus-col", default="locus",
                   help="--pop-freqs / --write-pop-freqs: locus column (default: locus)")
    p.add_argument("--freq-allele-col", default="allele",
                   help="--pop-freqs / --write-pop-freqs: allele column (default: allele)")
    p.add_argument("--freq-col", default="freq",
                   help="--pop-freqs / --write-pop-freqs: frequency column (default: freq; "
                        "renormalised per locus, so counts work too)")
    p.add_argument("--n-bins", type=int, default=10,
                   help="Number of MAF bins (default: 10, moimix). 0 = no binning: regress "
                        "every site's own Hw on its Hs. Use 0 for --allele-table, where "
                        "loci often have no allele above 50%% and the MAF bins do not fit.")
    p.add_argument("--min-alt-samples", type=int, default=0,
                   help="Keep only sites with the alt seen in >= this many samples "
                        "(default: 0; CNV gate uses 2)")
    p.add_argument("--no-snps-only", dest="snps_only", action="store_false",
                   help="Score on every record with AD (indels, MNPs), not only sites whose "
                        "REF and every ALT are single bases (the default)")
    p.add_argument("--snps-only", dest="snps_only", action="store_true", default=True,
                   help=argparse.SUPPRESS)  # the default; kept so older scripts still run
    p.add_argument("--multiallelic", choices=MULTIALLELIC_MODES, default="collapse",
                   help="What to do with a record that has more than one ALT: 'collapse' "
                        "(default) counts every allele's depth, so heterozygosity is "
                        "1 - sum(p^2) over all alleles and a mix of two ALTs reads as "
                        "mixed; 'skip' drops those records (what moimix::getFws sees on "
                        "a `bcftools norm -m-` split callset). At biallelic sites the two "
                        "are identical.")
    p.add_argument("--no-trim", dest="trim", action="store_false",
                   help="Keep ALT alleles that no sample here has reads for. By default they "
                        "are dropped before a site is classified (bcftools view "
                        "--trim-alt-alleles, done on AD), so a callset joint-called across "
                        "a larger cohort does not turn SNP sites into indel sites or "
                        "biallelic ones into multiallelic ones. Needed for moimix parity.")
    p.add_argument("--monoclonal-threshold", type=float, default=0.95,
                   help="Fws >= this is reported monoclonal (default: 0.95)")
    p.add_argument("--population-name",
                   help="If set, add a population_name column with this value (eases merging)")
    p.add_argument("--exclude-call-regions",
                   help="TSV (chrom, call_start, call_end) of CNV windows to exclude from Fws")
    p.add_argument("--out", default="-", help="Output TSV ('-' = STDOUT, default)")
    return p


def parse_args_calculate_fws():
    return get_parser_calculate_fws().parse_args()


def _load_samples(spec):
    sp = Path(spec)
    if sp.exists():
        return [x.strip() for x in sp.read_text().splitlines() if x.strip()]
    return [x for x in spec.split(",") if x]


def calculate_fws():
    args = parse_args_calculate_fws()
    exclude = load_exclude_regions(args.exclude_call_regions)

    if args.allele_table:
        if exclude:
            sys.exit("--exclude-call-regions does not apply to --allele-table (no positions)")
        samples, depths = read_allele_table(
            args.allele_table, sample_col=args.sample_col, locus_col=args.locus_col,
            allele_col=args.allele_col, reads_col=args.reads_col)
        if args.n_bins != 0:
            print(f"fws: note: --allele-table with --n-bins {args.n_bins}; microhaplotype "
                  f"loci often have no allele above 50%, so --n-bins 0 (per-locus "
                  f"regression) is the recommended setting", file=sys.stderr)
    elif args.input_vcf:
        samples, depths = read_ad_vcf(args.input_vcf, exclude, snps_only=args.snps_only,
                                      multiallelic=args.multiallelic, trim=args.trim)
    else:
        if not args.samples:
            sys.exit("--ad-table requires --samples (the column order)")
        samples = _load_samples(args.samples)
        depths = read_ad_table(args.ad_table, samples, exclude, snps_only=args.snps_only,
                               multiallelic=args.multiallelic, trim=args.trim)

    if depths.n_sites == 0:
        sys.exit("no usable loci found in the allele table" if args.allele_table else
                 "no usable " + ("SNP" if args.snps_only else "variant")
                 + " sites with AD found in the input")

    fcols = dict(locus_col=args.freq_locus_col, allele_col=args.freq_allele_col,
                 freq_col=args.freq_col)
    if args.write_pop_freqs:
        write_pop_freqs(depths.population_freqs(), args.write_pop_freqs, **fcols)
    pop_freqs = read_pop_freqs(args.pop_freqs, **fcols) if args.pop_freqs else None
    if pop_freqs is not None and not pop_freqs:
        sys.exit(f"--pop-freqs {args.pop_freqs}: no usable rows (check --freq-*-col)")

    fws, n_info = compute_fws(depths, estimator=args.estimator, min_depth=args.min_depth,
                              n_bins=args.n_bins, min_alt_samples=args.min_alt_samples,
                              pop_freqs=pop_freqs)
    if pop_freqs is not None:
        miss = compute_fws.last_pop_freq_misses
        if miss == depths.n_sites:
            sys.exit(f"--pop-freqs: none of the {depths.n_sites} sites in the input has an "
                     f"entry (loci are keyed CHROM:POS for VCF/AD input, by name for "
                     f"--allele-table)")
        print(f"fws: population frequencies from {args.pop_freqs}: {depths.n_sites - miss} "
              f"of {depths.n_sites} sites matched" + (f", {miss} dropped" if miss else ""),
              file=sys.stderr)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        header = ["sample", "fws", "n_sites", "monoclonal"]
        if args.population_name is not None:
            header.append("population_name")
        out.write("\t".join(header) + "\n")
        for s, f, n in zip(samples, fws, n_info):
            mono = "" if not np.isfinite(f) else str(bool(f >= args.monoclonal_threshold))
            fval = "" if not np.isfinite(f) else f"{f:.6f}"
            row = [s, fval, str(int(n)), mono]
            if args.population_name is not None:
                row.append(args.population_name)
            out.write("\t".join(row) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()

    n_ok = int(np.isfinite(fws).sum())
    n_mono = int(np.nansum(fws >= args.monoclonal_threshold))
    note = depths.multiallelic_note()
    print(f"fws: {n_ok} samples scored over {depths.n_sites} sites ({args.estimator}); "
          f"{n_mono} monoclonal (Fws>={args.monoclonal_threshold})"
          + (f"; {note}" if note else ""), file=sys.stderr)


if __name__ == "__main__":
    calculate_fws()
