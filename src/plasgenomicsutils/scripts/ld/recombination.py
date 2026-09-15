#!/usr/bin/env python
"""Population recombination-rate map from the monoclonal isolates, via LDhat or pyrho.

We do the part specific to a haploid malaria callset -- select the monoclonal isolates,
pull each out as a single haplotype, drop the sites the estimators cannot use -- and hand
the clean phased panel to a coalescent composite-likelihood estimator (``--emit``):

* ``ldhat``  -- **recommended for Pf.** LDhat's ``interval`` rjMCMC resolves gene-scale
  structure at these SNP densities. On Pf WGS, pyrho collapses each chromosome to ~one
  rate while LDhat recovers real cold-/hot-spots (e.g. pfpx1, pfkelch13 read as cold-spots).
* ``pyrho``  -- kept as a fast alternative, but its fused-LASSO over-smooths Pf data.

Why monoclonal only: both estimators read every haplotype as one phased sequence. A mixed
infection (COI>1, Fws below the bar) has no defined phase, so its "haplotype" is a blend of
lineages and the two-locus LD signal rho is estimated from is corrupted. The Fws gate is
the modelling assumption, not a nicety.

Neither estimator is imported -- both are external tools we shell out to. pyrho lives in
its own conda env (deps that do not co-resolve here); LDhat is a separate C build (see
``scripts/build_ldhat.sh``; on Apple Silicon it must be x86_64, run under Rosetta). See
``--pyrho-cmd`` / ``--ldhat-dir``. Note the rho conventions differ: pyrho reports 2*Ne*r
per bp, LDhat 4*Ne*r per bp (each output stamps its own in the TSV header).
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
from pathlib import Path

import numpy as np

from ...lib.ld import read_dosages
from ...lib.fws import monoclonal_samples
from ...utils.small_utils import Utils

DEFAULT_FWS_MIN = 0.95        # same bar as fws_filter / calculate_fws --monoclonal-threshold
DEFAULT_MISSING_MAX = 0.10    # per-site missing-rate ceiling before a site is dropped


# --------------------------------------------------------------------------- interface
def get_parser_ld_recombination() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasgenomicsutils ld_recombination",
        description="Per-SNP population recombination-rate map along each chromosome from "
                    "the monoclonal isolates, via LDhat (--emit ldhat, recommended for Pf) "
                    "or pyrho (--emit pyrho). Only monoclonal samples are used: each "
                    "haplotype is read as one phased sequence, and a mixed infection has no "
                    "phase to read.",
        epilog="Prerequisites depend on --emit. ldhat: the interval/lkgen/stat binaries "
               "(--ldhat-dir or $LDHAT_DIR; build with scripts/build_ldhat.sh). pyrho: a "
               "lookup table for a demography and a sample size >= the haplotype count "
               "(pyrho make_table --approx; hyperparam to pick --block-penalty/--window-"
               "size), run from its own env (--pyrho-cmd). Multiallelic and spanning-"
               "deletion sites are handled exactly as ld_decay handles them.\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # --- input / sample selection ---
    p.add_argument("--vcf", required=True,
                   help="Biallelic-SNP VCF/BCF (the full callset; do not LD-prune)")
    p.add_argument("--fws-table", default=None,
                   help="calculate_fws output (sample, fws, ...). Without it, --samples or "
                        "all samples are used and you assert they are monoclonal")
    p.add_argument("--fws-min", type=float, default=DEFAULT_FWS_MIN,
                   help="Keep isolates with Fws >= this (default: %(default)s)")
    p.add_argument("--samples", default=None,
                   help="Explicit monoclonal sample list (one per line); overrides the gate")
    p.add_argument("--region", action="append", default=None, metavar="REG",
                   help="chrom or chrom:start-end; repeatable (needs an index)")
    p.add_argument("--min-depth", type=int, default=5,
                   help="Treat a call backed by fewer reads as missing (default: %(default)s)")
    p.add_argument("--max-missing", type=float, default=DEFAULT_MISSING_MAX,
                   help="Drop a site whose missing rate exceeds this; pyrho needs complete "
                        "columns, and the rest are imputed to the major allele "
                        "(default: %(default)s)")
    p.add_argument("--maf", type=float, default=0.0,
                   help="Minor-allele-frequency floor among the monoclonal panel "
                        "(default: %(default)s, keep all)")
    # --- pyrho ---
    p.add_argument("--table", default=None,
                   help="pyrho lookup table (.hdf) for a sample size >= the haplotype "
                        "count. Required for --emit pyrho")
    p.add_argument("--block-penalty", type=float, default=50.0,
                   help="pyrho fused-LASSO block penalty (from pyrho hyperparam; "
                        "default: %(default)s)")
    p.add_argument("--window-size", type=int, default=50,
                   help="pyrho window size in SNPs (from pyrho hyperparam; "
                        "default: %(default)s)")
    p.add_argument("--numthreads", type=int, default=1)
    p.add_argument("--pyrho-cmd", default="auto",
                   help="How to invoke pyrho. 'auto' (default): use `pyrho` on PATH, else "
                        "`mamba run -n pyrho pyrho`. Or give an explicit command, e.g. "
                        "'mamba run -n myenv pyrho' or a full path")
    # --- LDhat (the recommended method for Pf: pyrho's fused-LASSO over-smooths this
    # SNP density -- on Pf WGS it collapses each chromosome to ~1 rate, while LDhat's
    # interval rjMCMC resolves gene-scale structure at the same n) ---
    p.add_argument("--ldhat-dir", default=os.environ.get("LDHAT_DIR"),
                   help="Directory holding the LDhat `interval`, `lkgen` and `stat` "
                        "binaries (or set $LDHAT_DIR). Required for --emit ldhat")
    p.add_argument("--ldhat-lk", default=None,
                   help="Base LDhat likelihood lookup table (plain or .gz) to subsample to "
                        "the panel size with lkgen. Default: <ldhat-dir>/lk_files/"
                        "lk_n120_t0.001.gz. Its sequence count must be >= the haplotype count")
    p.add_argument("--ldhat-arch-prefix", default="auto",
                   help="Command prefix to run the LDhat binaries under. 'auto' uses "
                        "`arch -x86_64` on Apple Silicon (the LDhat C bus-errors when built "
                        "native arm64) and nothing elsewhere; 'none' forces no prefix")
    p.add_argument("--ldhat-bpen", type=float, default=5.0,
                   help="LDhat interval block penalty (default: %(default)s, as Niare et al.)")
    p.add_argument("--ldhat-its", type=int, default=5_000_000,
                   help="LDhat interval MCMC iterations (default: %(default)s)")
    p.add_argument("--ldhat-samp", type=int, default=5000,
                   help="LDhat interval sampling interval (default: %(default)s)")
    p.add_argument("--ldhat-burn", type=float, default=0.2,
                   help="Fraction of LDhat samples discarded as burn-in (default: %(default)s)")
    p.add_argument("--max-snps-per-block", type=int, default=3000,
                   help="Split a chromosome above this many SNPs into overlapping LDhat "
                        "blocks and stitch (interval's rjMCMC cost grows with SNP count; "
                        "0 disables splitting) (default: %(default)s)")
    p.add_argument("--block-overlap", type=int, default=200,
                   help="SNPs of overlap between adjacent LDhat blocks; half is trimmed "
                        "from each interior edge when stitching (default: %(default)s)")
    # --- outputs / control ---
    p.add_argument("--emit", choices=("pyrho", "ldhat", "ldhat-inputs", "haps"),
                   default="pyrho",
                   help="pyrho: run pyrho and write the rate map. ldhat: run LDhat "
                        "(interval + stat, block-stitched) and write the rate map -- the "
                        "recommended method for Pf. ldhat-inputs: write LDhat sites/locs and "
                        "stop. haps: write the 0/1 panel + positions and stop -- no external "
                        "tool needed (default: %(default)s)")
    p.add_argument("--out-prefix", required=True,
                   help="Prefix for the rate-map TSV and any intermediate input files")
    p.add_argument("--keep-intermediate", action="store_true",
                   help="Keep the per-chromosome intermediates (pyrho input VCFs, or the "
                        "LDhat block workdir) instead of deleting them")
    p.add_argument("--overwrite", action="store_true")
    return p


def parse_args_ld_recombination():
    return get_parser_ld_recombination().parse_args()


# ----------------------------------------------------------------- pyrho command resolve
def resolve_pyrho_cmd(spec):
    """Return the argv prefix that runs pyrho. 'auto' prefers `pyrho` on PATH, then falls
    back to `mamba run -n pyrho pyrho` (or conda). An explicit spec is split on spaces."""
    if spec and spec != "auto":
        return spec.split()
    if shutil.which("pyrho"):
        return ["pyrho"]
    for mgr in ("mamba", "micromamba", "conda"):
        if shutil.which(mgr):
            return [mgr, "run", "-n", "pyrho", "pyrho"]
    raise SystemExit(
        "pyrho not found: put it on PATH, build the 'pyrho' env (see environment.yml), "
        "or pass --pyrho-cmd explicitly")


# --------------------------------------------------------------- the pyrho input step
def build_haplotype_panel(gn, chrom, pos, *, max_missing, maf):
    """Turn the (n_variants, n_samples) 0/2/-1 dosage array from ``read_dosages`` into a
    complete, biallelic 0/1 haplotype panel pyrho / LDhat can consume.

    ``read_dosages`` already dropped multiallelics and masked spanning-deletion calls to
    -1, coding hom-ref 0 / hom-alt 2 / missing -1 (het is -1 under het="missing"). For a
    haploid that 0/2 *is* the haplotype -- we rescale to 0/1, enforce completeness, and
    yield one chromosome at a time.

    Yields ``(chrom, hap, positions)`` where ``hap`` is int8 ``(n_snps, n_haplotypes)``
    0/1, no missing, no monomorphic sites, and ``positions`` are 0-based.
    """
    hap01 = np.where(gn == 2, 1, gn).astype(np.int8)   # 0->0, 2->1, -1 stays -1
    for c in _unique_stable(chrom):
        m = chrom == c
        sub, p = hap01[m], pos[m]
        miss = sub < 0
        keep = miss.mean(axis=1) <= max_missing
        sub, p, miss = sub[keep], p[keep], miss[keep]
        if miss.any():                    # impute the rest to the site major allele
            major = (np.where(miss, 0, sub).sum(axis=1) >=
                     (~miss).sum(axis=1) / 2.0).astype(np.int8)
            sub = np.where(miss, major[:, None], sub)
        alt_freq = sub.mean(axis=1)
        maf_site = np.minimum(alt_freq, 1 - alt_freq)
        poly = (maf_site > 0) & (maf_site >= maf)     # rho undefined at a fixed site
        if poly.any():
            yield c, sub[poly], p[poly]


def _unique_stable(a):
    seen, out = set(), []
    for x in a:
        if x not in seen:
            seen.add(x); out.append(x)
    return out


def write_pyrho_vcf(chrom, hap, positions, path):
    """Write the phased VCF pyrho's reader wants.

    pyrho requires *diploid* GT columns and, under --ploidy 1, flattens each ``a|b`` back
    into two independent haplotypes -- so we pair consecutive haplotypes into phased
    pseudo-diploids. The pairing is arbitrary and lossless: flattening reproduces exactly
    the original haplotype set. An odd panel drops its last haplotype (pyrho needs pairs).

    Returns the number of haplotypes actually written.
    """
    n_hap = hap.shape[1]
    dropped_odd = n_hap % 2
    if dropped_odd:
        hap, n_hap = hap[:, :-1], n_hap - 1
    n_pairs = n_hap // 2
    cols = [f"pd{i}" for i in range(n_pairs)]
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={chrom}>\n")         # silences pyrho's "contig not defined"
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(cols) + "\n")
        for j in range(hap.shape[0]):
            row = hap[j]
            gt = "\t".join(f"{int(row[2*k])}|{int(row[2*k+1])}" for k in range(n_pairs))
            fh.write(f"{chrom}\t{positions[j] + 1}\t.\tA\tT\t.\tPASS\t.\tGT\t{gt}\n")
    return n_hap, dropped_odd


def write_ldhat(chrom, hap, positions, sites_path, locs_path):
    """Classic LDhat inputs, for anyone who wants canonical LDhat numbers instead.

    LDhat's ``interval``/``pairwise`` require STRICTLY increasing positions in the locs
    file (kb). Two records at the same site -- a SNP and a co-located record the panel did
    not merge -- or a pair that collides at the 3-decimal (1 bp) kb resolution otherwise
    trips "SNPs must be monotonically increasing". We nudge any tie up by 1 bp so the file
    is always accepted; the displacement is below LDhat's own resolution.
    """
    n_snp, n_hap = hap.shape
    with open(sites_path, "w") as fh:
        fh.write(f"{n_hap} {n_snp} 1\n")            # haploid -> 1
        for i in range(n_hap):
            fh.write(f">hap{i}\n" + "".join(str(int(x)) for x in hap[:, i]) + "\n")
    kb = (positions / 1000.0).astype(float).copy()
    for j in range(1, len(kb)):                    # enforce strict monotonicity
        if kb[j] <= kb[j - 1]:
            kb[j] = kb[j - 1] + 0.001              # +1 bp
    with open(locs_path, "w") as fh:
        fh.write(f"{n_snp} {kb[-1]:.3f} L\n" + " ".join(f"{x:.3f}" for x in kb) + "\n")


def run_pyrho_optimize(vcf_path, out_path, *, cmd, table, block_penalty, window_size,
                       numthreads):
    argv = list(cmd) + ["optimize", "--vcffile", str(vcf_path), "--tablefile", str(table),
                        "--ploidy", "1", "--blockpenalty", str(block_penalty),
                        "--windowsize", str(window_size), "--numthreads", str(numthreads),
                        "--outfile", str(out_path)]
    print("  " + " ".join(argv))
    subprocess.run(argv, check=True)


def read_pyrho_map(path, chrom):
    """pyrho optimize writes headerless ``start<TAB>end<TAB>rho`` rows (rho per bp between
    start and end). Return a DataFrame with the chromosome prepended."""
    import pandas as pd
    df = pd.read_csv(path, sep="\t", header=None, names=["start", "end", "rho_per_bp"])
    df.insert(0, "chrom", chrom)
    return df


# ------------------------------------------------------------------ LDhat interval runner
def resolve_ldhat(ldhat_dir, arch_prefix="auto"):
    """Return ``(prefix, {name: path})`` for the LDhat binaries. ``prefix`` is prepended to
    every invocation -- ``arch -x86_64`` on Apple Silicon, where the LDhat C bus-errors when
    built native arm64 and must run under Rosetta."""
    if not ldhat_dir:
        raise SystemExit("--emit ldhat needs --ldhat-dir (or $LDHAT_DIR): the folder with "
                         "the interval/lkgen/stat binaries")
    d = Path(ldhat_dir)
    bins = {}
    for name in ("interval", "lkgen", "stat"):
        b = d / name
        if not b.exists():
            raise SystemExit(f"LDhat '{name}' binary not found in {ldhat_dir}")
        bins[name] = str(b)
    if arch_prefix == "none":
        prefix = []
    elif arch_prefix == "auto":
        prefix = (["arch", "-x86_64"] if platform.system() == "Darwin"
                  and platform.machine() == "arm64" else [])
    else:
        prefix = arch_prefix.split()
    return prefix, bins


def _lookup_nseq(path):
    """First whitespace token of an LDhat lookup file is its sequence count."""
    opener = __import__("gzip").open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line:
                return int(line.split()[0])
    raise SystemExit(f"{path}: empty LDhat lookup file")


def prep_ldhat_lookup(base_lk, n, workdir, *, prefix, lkgen):
    """Subsample the base LDhat lookup to ``n`` sequences with lkgen, cached per n. lkgen
    reads the table as plain text, so a ``.gz`` base is gunzipped first."""
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True)
    cached = workdir / f"new_lk_n{n}.txt"
    if cached.exists() and cached.stat().st_size > 0:
        return str(cached)
    base = Path(base_lk)
    if str(base).endswith(".gz"):
        plain = workdir / base.with_suffix("").name
        if not plain.exists() or plain.stat().st_size == 0:
            import gzip
            with gzip.open(base, "rb") as fi, open(plain, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        base = plain
    base_n = _lookup_nseq(base)
    if base_n < n:
        raise SystemExit(f"LDhat lookup {base} has n={base_n} sequences but the panel has "
                         f"{n}; supply a larger --ldhat-lk (lkgen only downsamples)")
    subprocess.run(prefix + [lkgen, "-lk", str(base), "-nseq", str(n)],
                   cwd=str(workdir), check=True, stdout=subprocess.DEVNULL)
    gen = workdir / "new_lk.txt"
    if not gen.exists() or gen.stat().st_size == 0:
        raise SystemExit("lkgen produced no usable new_lk.txt")
    if _lookup_nseq(gen) != n:
        raise SystemExit("lkgen output has the wrong sequence count")
    gen.rename(cached)
    return str(cached)


def parse_ldhat_res(path):
    """LDhat ``stat`` res.txt -> ``(kb, rho)``. Columns are locus position (kb) and the
    per-interval rate Mean_rho = 4*Ne*r per kb (the first row, locus -1, is the map total
    and is dropped)."""
    kb, rho = [], []
    for line in Path(path).read_text().splitlines()[1:]:
        f = line.split()
        if len(f) >= 2 and float(f[0]) >= 0:
            kb.append(float(f[0])); rho.append(float(f[1]))
    return np.array(kb), np.array(rho)


def run_ldhat_block(chrom, hap, positions, *, lookup, prefix, bins, bpen, its, samp,
                    burn, seed, workdir):
    """Run interval + stat on one panel and return ``(pos_bp, rho_per_kb)`` per locus."""
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True)
    sites, locs = str(workdir / "sites.txt"), str(workdir / "locs.txt")
    write_ldhat(chrom, hap, positions, sites, locs)
    subprocess.run(prefix + [bins["interval"], "-seq", sites, "-loc", locs, "-lk", lookup,
                             "-its", str(its), "-bpen", str(bpen), "-samp", str(samp),
                             "-seed", str(seed)],
                   cwd=str(workdir), check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    n_burn = int(burn * its / samp)
    subprocess.run(prefix + [bins["stat"], "-input", "rates.txt", "-loc", locs,
                             "-burn", str(n_burn)],
                   cwd=str(workdir), check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    kb, rho = parse_ldhat_res(workdir / "res.txt")
    return kb * 1000.0, rho          # kb positions -> bp (absolute, since locs used abs pos)


def run_ldhat(chrom, hap, positions, *, lookup, prefix, bins, bpen, its, samp, burn,
              max_snps, overlap, workdir, seed=7):
    """Estimate a recombination map for one chromosome, splitting into overlapping blocks
    above ``max_snps`` SNPs and stitching (interval's rjMCMC cost grows with SNP count).

    Returns a DataFrame ``chrom, start, end, rho_per_bp`` where rho is 4*Ne*r per bp.
    """
    import pandas as pd
    n_snp = hap.shape[0]
    step = max(1, max_snps - overlap)
    blocks = ([(0, n_snp)] if (max_snps <= 0 or n_snp <= max_snps)
              else [(s, min(s + max_snps, n_snp)) for s in range(0, n_snp, step)
                    if s < n_snp])
    pos_all, rho_all = [], []
    for bi, (s, e) in enumerate(blocks):
        bdir = Path(workdir) / f"block_{bi:02d}"
        bp, rho = run_ldhat_block(chrom, hap[s:e], positions[s:e], lookup=lookup,
                                  prefix=prefix, bins=bins, bpen=bpen, its=its, samp=samp,
                                  burn=burn, seed=seed, workdir=bdir)
        # trim half the overlap from interior edges so blocks tile without double-counting
        lo = positions[s + overlap // 2] if bi > 0 else -np.inf
        hi = positions[e - 1 - overlap // 2] if bi < len(blocks) - 1 else np.inf
        keep = (bp >= lo) & (bp <= hi)
        pos_all.append(bp[keep]); rho_all.append(rho[keep])
    pos = np.concatenate(pos_all); rho = np.concatenate(rho_all)
    order = np.argsort(pos, kind="stable")
    pos, rho = pos[order], rho[order]
    uniq = np.concatenate([[True], np.diff(pos) > 0])       # drop any overlap duplicates
    pos, rho = pos[uniq], rho[uniq]
    # each locus carries the rate of the interval to its right
    start = pos[:-1].astype(np.int64)
    end = pos[1:].astype(np.int64)
    return pd.DataFrame({"chrom": chrom, "start": start, "end": end,
                         "rho_per_bp": rho[:-1] / 1000.0})


# ---------------------------------------------------------------------------- driver
def ld_recombination():
    args = parse_args_ld_recombination()
    rate_map = f"{args.out_prefix}.recomb_map.tsv"
    pyrho_cmd = ldhat = None
    if args.emit == "pyrho":
        if not args.table:
            raise SystemExit("--emit pyrho needs --table (a pyrho make_table .hdf)")
        Utils.output_file_check(rate_map, args.overwrite)
        pyrho_cmd = resolve_pyrho_cmd(args.pyrho_cmd)
    elif args.emit == "ldhat":
        Utils.output_file_check(rate_map, args.overwrite)
        ldhat_prefix, ldhat_bins = resolve_ldhat(args.ldhat_dir, args.ldhat_arch_prefix)
        base_lk = args.ldhat_lk or str(Path(args.ldhat_dir) / "lk_files"
                                       / "lk_n120_t0.001.gz")
        if not Path(base_lk).exists():
            raise SystemExit(f"LDhat lookup table not found: {base_lk} (pass --ldhat-lk)")
        ldhat_work = Path(f"{args.out_prefix}.ldhat_work")
        ldhat = (ldhat_prefix, ldhat_bins, base_lk, ldhat_work)

    if args.samples:
        keep = [s.strip() for s in Path(args.samples).read_text().split() if s.strip()]
        print(f"[gate] {len(keep)} sample(s) from --samples")
    elif args.fws_table:
        keep = monoclonal_samples(args.fws_table, fws_min=args.fws_min)
        print(f"[gate] {len(keep)} isolate(s) with Fws >= {args.fws_min:g} in the table")
    else:
        keep = None
        print("[warn] no --fws-table/--samples: using every sample and ASSUMING monoclonal")

    gn, chrom, pos, names, _counts = read_dosages(
        args.vcf, samples=keep, regions=args.region,
        het="missing", min_depth=args.min_depth)
    print(f"[read] {gn.shape[0]:,} biallelic SNPs x {gn.shape[1]} monoclonal isolate(s)")

    per_chrom = list(build_haplotype_panel(
        gn, chrom, pos, max_missing=args.max_missing, maf=args.maf))
    if not per_chrom:
        raise SystemExit("no site survived completeness/MAF filtering")

    ldhat_lookup = None
    frames = []
    for c, hap, p in per_chrom:
        stem = f"{args.out_prefix}.{c}"
        if args.emit == "haps":
            np.savetxt(f"{stem}.hap01.tsv", hap.T, fmt="%d")   # rows = haplotypes
            np.savetxt(f"{stem}.pos.tsv", p, fmt="%d")
            print(f"  [{c}] {hap.shape[0]:,} SNPs x {hap.shape[1]} haplotypes -> {stem}.hap01.tsv")
            continue
        if args.emit == "ldhat-inputs":
            write_ldhat(c, hap, p, f"{stem}.sites", f"{stem}.locs")
            print(f"  [{c}] wrote {stem}.sites / {stem}.locs")
            continue
        if args.emit == "ldhat":
            prefix, bins, base_lk, work = ldhat
            if ldhat_lookup is None:                 # subsample the lookup once, cached
                ldhat_lookup = prep_ldhat_lookup(base_lk, hap.shape[1], work,
                                                 prefix=prefix, lkgen=bins["lkgen"])
            n_snp = hap.shape[0]
            n_blk = 1 if (args.max_snps_per_block <= 0 or n_snp <= args.max_snps_per_block) \
                else -(-n_snp // max(1, args.max_snps_per_block - args.block_overlap))
            print(f"  [{c}] LDhat interval on {n_snp:,} SNPs x {hap.shape[1]} hap"
                  f"{'' if n_blk == 1 else f' in {n_blk} blocks'} (bpen {args.ldhat_bpen:g})")
            frames.append(run_ldhat(c, hap, p, lookup=ldhat_lookup, prefix=prefix, bins=bins,
                                    bpen=args.ldhat_bpen, its=args.ldhat_its,
                                    samp=args.ldhat_samp, burn=args.ldhat_burn,
                                    max_snps=args.max_snps_per_block,
                                    overlap=args.block_overlap, workdir=work / c))
            continue
        vcf_in = f"{stem}.pyrho_in.vcf"
        n_hap, odd = write_pyrho_vcf(c, hap, p, vcf_in)
        if odd:
            print(f"  [{c}] odd panel: dropped 1 haplotype for pairing (n={n_hap})")
        out = f"{stem}.pyrho_out.tsv"
        run_pyrho_optimize(vcf_in, out, cmd=pyrho_cmd, table=args.table,
                           block_penalty=args.block_penalty, window_size=args.window_size,
                           numthreads=args.numthreads)
        frames.append(read_pyrho_map(out, c))
        if not args.keep_intermediate:
            Path(vcf_in).unlink(missing_ok=True)

    if args.emit in ("pyrho", "ldhat"):
        import pandas as pd
        full = pd.concat(frames, ignore_index=True)
        if args.emit == "pyrho":
            header = (f"#tool=pyrho\t#ploidy=1\t#fws_min={args.fws_min}\t"
                      f"#block_penalty={args.block_penalty}\t#window_size={args.window_size}\t"
                      f"#note=rho is 2*Ne*r per bp for a haploid; divide by 2*Ne for r/bp")
        else:
            header = (f"#tool=ldhat\t#fws_min={args.fws_min}\t#bpen={args.ldhat_bpen}\t"
                      f"#its={args.ldhat_its}\t#maf={args.maf}\t"
                      f"#note=rho is 4*Ne*r per bp (LDhat convention); divide by 4*Ne for r/bp")
        Utils.write_tsv_gz(full, rate_map, header_comment=header)
        print(f"[done] wrote {rate_map} ({len(full)} intervals across {len(frames)} chrom)")
        if args.emit == "ldhat" and not args.keep_intermediate:
            shutil.rmtree(ldhat[3], ignore_errors=True)
    else:
        print(f"[done] wrote {args.emit} inputs with prefix {args.out_prefix}")


if __name__ == "__main__":
    ld_recombination()
