"""EIGENSTRAT-style IBD selection statistic (XiR,s).

Allele frequencies must be supplied externally and are never proxied from the
IBD matrix.

Method:
  1. binary IBD matrix (pairs x SNPs)
  2. subtract per-pair means -- removes each pair's overall relatedness
  3. divide by sqrt(p(1-p)), p = SNP allele frequency
  4. sum over pairs / sqrt(n_pairs) -> raw per-SNP statistic
  5. bin SNPs into equal-frequency MAF bins; within-bin z-score
  6. z -> p -> -log10(p), upper tail by default (see ``tail``)
  7. multiple-testing views, and the diagnostics that say what each is worth

**Two variants, and why.** Henden et al. (PLoS Genet 2018) describe XiR,s with a second
centring between steps 2 and 4: *"we subtract the row mean from each row"*, each row being
one SNP, followed by *"we calculate row sums"*. Centring a row and then summing that same
row cancels exactly, so the described statistic is identically zero and what survives is
floating-point residue. That residue is not random -- it grows with the number of pairs
sharing, so it behaves like a noisy, uncalibrated proxy for excess sharing, which is why
the recipe has looked serviceable. It is not reproducible: at float32 the residue is ~1e-5,
at float64 ~1e-14, and the two rank SNPs almost independently. The same cancellation is
present in isoRelate's ``iRfunction`` and in ibdutils' ``calc_xirs_raw_stats_per_chr``; on
one real cohort those two and this one agreed on 29 of their top 100 SNPs.

``variant="corrected"`` (the default) omits the per-SNP centring, so step 4 measures what
the method is for: sharing at a locus above what these pairs' relatedness predicts.
``variant="published"`` restores the cancellation, in float32, purely to reproduce output
from earlier versions of this tool; it prints a warning and cannot be made to agree with
isoRelate or ibdutils, which cancel at a different precision.

**Tail.** ``tail="upper"`` (the default) asks only whether a locus is shared *more* than
expected, which is what a positive-selection scan claims. ``tail="two-sided"`` squares the
z-score into a chi2(1), which scores a sharing *deficit* exactly like an excess -- what the
published recipe does. Every row carries a ``direction`` column either way.

**Choosing a threshold.** Three are reported, and they differ in what they control and in
what they assume:

* **Bonferroni** (`significant`) -- family-wise: near-certainty that no SNP called is a
  false positive. Reads the p-value off a chi2(1).
* **Benjamini-Hochberg** (`significant_fdr`) -- a stated share of the SNPs called may be
  false, so it calls more. Also off the chi2(1).
* **Permutation** (`significant_perm`, `significant_fdr_perm`) -- from
  :func:`permutation_null`, which draws the reference from the data instead. Family-wise
  and FDR flavours both.

`lambda_gc` decides which to trust: the median chi2 over the chi2(1) median, 1 when the
reference fits. The z-scores are standardised to zero mean and unit variance *within each
MAF bin*, which fixes the first two moments and leaves the shape alone -- and the shape of
IBD sharing is nothing like a normal. On real *P. falciparum* data lambda comes out near
0.1, a tight bulk with very heavy tails, so the two chi2(1) thresholds are miscalibrated
and the permutation pair is the one to quote.

Either way the *ranking* stands: every step from `z_score` to `neg_log10_p` is monotone, so
a poorly fitting reference mislabels the axis without moving any SNP relative to another.
`pval` and `q_value` inherit the misfit and should not be read as probabilities;
`p_empirical` / `q_empirical` are their calibrated counterparts.

No correction here knows about linkage. Adjacent SNPs in one sweep are not independent
tests, so SNP counts overstate the number of findings; merge significant SNPs into peaks
before counting discoveries.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

from .intervals import SNP_COORD_SYSTEM, check_snp_coord_system
from ..utils.small_utils import Utils


# ---------------------------------------------------------------------------
# label / AF loading
# ---------------------------------------------------------------------------

def parse_snp_labels(snp_labels: list, with_pos_vcf: bool = False) -> pd.DataFrame:
    """Split ``chr:pos0`` labels into columns. ``pos`` is 0-based, like the label.

    ``with_pos_vcf`` adds the 1-based VCF position as ``pos_vcf`` for cross-referencing a
    variant against the VCF or a browser; it is redundant with ``pos``, so it is off by
    default rather than inflating every per-SNP table.
    """
    rows = []
    for label in snp_labels:
        if ":" in label:
            chrom, pos = label.rsplit(":", 1)
            rows.append({"snp_id": label, "chr": chrom, "pos": int(pos)})
        else:
            rows.append({"snp_id": label, "chr": "unknown", "pos": -1})
    df = pd.DataFrame(rows)
    if with_pos_vcf and len(df):
        df["pos_vcf"] = np.where(df["pos"] >= 0, df["pos"] + 1, -1)
    return df


def _read_af_table(af_path: str, usecols: list) -> pd.DataFrame:
    """Read an AF table, verifying its stamped SNP coordinate system first."""
    with Utils.smart_open_read(af_path) as fh:
        first = fh.readline().strip()
    stamped = first.split("=", 1)[1] if first.startswith("#snp_coord_system=") else None
    check_snp_coord_system(stamped, af_path)
    try:
        return pd.read_csv(af_path, sep="\t", usecols=usecols, comment="#")
    except ValueError as e:
        have = pd.read_csv(af_path, sep="\t", nrows=0, comment="#").columns.tolist()
        missing = [c for c in usecols if c not in have]
        if not missing:
            raise
        raise SystemExit(
            f"ERROR: {af_path} has no column(s): {', '.join(missing)}.\n"
            f"  columns present: {', '.join(have)}\n"
            "  Tables written before this version hold only snp_id and af; re-run "
            "compute_allele_freqs to get the others."
        ) from e


def _reject_per_alt_table(df: pd.DataFrame, af_path: str, keys: list[str]) -> None:
    """Refuse a table with more than one row per SNP.

    ``compute_allele_freqs --per-alt`` writes one row per (SNP, ALT), and the two modes
    write the same filename. The joins below are dictionary lookups, so a per-ALT table
    would not fail -- it would silently keep whichever alternate happened to be last and
    run the whole statistic on that one allele's frequency. At a biallelic site the two
    modes agree row for row, which is why this only bites once the callset carries a
    multiallelic record.
    """
    dup = df.duplicated(subset=keys, keep=False)
    if not dup.any():
        return
    offenders = df.loc[dup, keys].drop_duplicates()
    example = " / ".join(f"{k}={v}" for k, v in offenders.iloc[0].items())
    raise SystemExit(
        f"ERROR: {af_path} has {len(offenders):,} {'/'.join(keys)} value(s) on more than "
        f"one row (e.g. {example}).\n"
        "  That is the shape `compute_allele_freqs --per-alt` writes: one row per "
        "(SNP, ALT).\n"
        "  The selection statistic joins one row per SNP, so re-run "
        "compute_allele_freqs without --per-alt."
    )


def load_global_af(af_path: str, snp_labels: list, af_col: str = "af") -> np.ndarray:
    """Global AFs aligned to ``snp_labels``; any missing SNP is a hard error.

    ``af_col`` names the column to read. The default is the allele-count frequency, which
    is the one that matches how the IBD matrix was built: hmmibd-rs reduces each sample to
    a single dominant allele, so the frequency of those same hard calls is what the
    expected-sharing model is about. ``af_weighted`` describes the within-host composition
    instead -- it counts minor clones the IBD analysis never saw -- so it answers a
    different question, and on the monoclonal cohorts this is usually run on the two agree
    anyway. Exposed so that can be checked rather than assumed.
    """
    af_df = _read_af_table(af_path, ["snp_id", af_col])
    _reject_per_alt_table(af_df, af_path, ["snp_id"])
    af_map = af_df.set_index("snp_id")[af_col].to_dict()
    af = np.array([af_map.get(s, np.nan) for s in snp_labels])
    missing = int(np.isnan(af).sum())
    if missing > 0:
        print(f"  ERROR: {missing:,} SNPs are missing from the AF file.")
        print("  First 10 SNP IDs expected (from matrix label file):")
        for s in snp_labels[:10]:
            print(f"    '{s}'")
        print("  First 10 SNP IDs found in AF file:")
        for s in list(af_map.keys())[:10]:
            print(f"    '{s}'")
        raise SystemExit(1)
    return af


def load_global_he(af_path: str, snp_labels: list) -> np.ndarray | None:
    """Expected heterozygosity aligned to ``snp_labels``, or ``None`` when the table has none.

    ``he`` is what generalises this statistic to any number of alleles, so it is read whenever
    the table carries it and the biallelic forms are used when it does not -- which is what a
    table written before the column existed requires. Returning ``None`` rather than raising
    keeps an old table working; the caller says which form it ended up using.
    """
    have = pd.read_csv(af_path, sep="\t", nrows=0, comment="#").columns.tolist()
    if "he" not in have:
        return None
    df = _read_af_table(af_path, ["snp_id", "he"])
    _reject_per_alt_table(df, af_path, ["snp_id"])
    m = df.set_index("snp_id")["he"].to_dict()
    return np.array([m.get(s, np.nan) for s in snp_labels])


def load_global_bin_inputs(af_path: str, snp_labels: list):
    """``(n_alleles_obs, maf_k)`` aligned to ``snp_labels``, or ``(None, None)``.

    The binning half of the k-allele generalisation, kept separate from ``he`` because they
    are different quantities used for different things -- see :func:`selection_bin_key`.
    Missing columns mean an older table, and the binning is then left exactly as it was.
    """
    have = pd.read_csv(af_path, sep="\t", nrows=0, comment="#").columns.tolist()
    if "maf_k" not in have or "n_alleles_obs" not in have:
        return None, None
    df = _read_af_table(af_path, ["snp_id", "n_alleles_obs", "maf_k"])
    _reject_per_alt_table(df, af_path, ["snp_id"])
    n = df.set_index("snp_id")["n_alleles_obs"].to_dict()
    k = df.set_index("snp_id")["maf_k"].to_dict()
    return (np.array([n.get(s, 2) for s in snp_labels]),
            np.array([k.get(s, np.nan) for s in snp_labels]))


def get_he_for_group(group, snp_labels, group_af_table, global_he):
    """Group ``he`` with the same fallback :func:`get_af_for_group` uses."""
    if (group_af_table is not None and "he" in group_af_table.columns
            and group in group_af_table["group"].values):
        sub = group_af_table[group_af_table["group"] == group]
        m = sub.set_index("snp_id")["he"].to_dict()
        return np.array([m.get(s, np.nan) for s in snp_labels])
    return None if global_he is None else global_he.copy()


def load_group_af_table(af_group_path: str, af_col: str = "af") -> pd.DataFrame:
    """Read the per-group allele-frequency table, verifying its coordinate stamp.

    ``af_col`` chooses which frequency column to use and arrives renamed to ``af``, so
    everything downstream joins on one name whichever was asked for.
    """
    want = ["group", "snp_id", af_col]
    have = pd.read_csv(af_group_path, sep="\t", nrows=0, comment="#").columns.tolist()
    if "he" in have and "he" not in want:
        want.append("he")
    df = _read_af_table(af_group_path, want)
    _reject_per_alt_table(df, af_group_path, ["group", "snp_id"])
    if af_col != "af":
        df = df.rename(columns={af_col: "af"})     # downstream joins on `af`
    # `he` rides along when the table has it -- it is what generalises the statistic to any
    # number of alleles -- but `af` stays the second-to-last column so the shape reads the
    # same whichever frequency column was asked for
    df = df[[c for c in ("group", "snp_id", "af") if c in df.columns]
            + [c for c in df.columns if c not in ("group", "snp_id", "af")]]
    print(f"  Loaded group AF table: {len(df):,} rows, "
          f"{df['group'].nunique()} groups, {df['snp_id'].nunique():,} SNPs")
    return df


def get_af_for_group(group, snp_labels, group_af_table, global_af) -> np.ndarray:
    """Group AF with fallback: per-group table -> global AF. Never from the matrix."""
    if group_af_table is not None and group in group_af_table["group"].values:
        sub = group_af_table[group_af_table["group"] == group]
        af_map = sub.set_index("snp_id")["af"].to_dict()
        af = np.array([af_map.get(s, np.nan) for s in snp_labels])
        missing = int(np.isnan(af).sum())
        if missing > 0:
            print(f"  ERROR: {missing:,} SNPs (out of {len(snp_labels):,}) missing from "
                  f"group AF file for group '{group}'.")
            raise SystemExit(1)
        return af
    return global_af.copy()


# ---------------------------------------------------------------------------
# pair subsetting
# ---------------------------------------------------------------------------

def within_group_row_indices(pair_labels, meta, group_col, group) -> np.ndarray:
    """Row indices for pairs where BOTH samples are from ``group``."""
    sample_to_group = meta.set_index("sample")[group_col].to_dict()
    idx = [
        i for i, label in enumerate(pair_labels)
        if all(sample_to_group.get(s) == group for s in label.split("__", 1))
    ]
    return np.array(idx, dtype=np.int64)


# ---------------------------------------------------------------------------
# core statistic
# ---------------------------------------------------------------------------

VARIANTS = ("corrected", "published")
TAILS = ("upper", "two-sided")

# The corrected statistic accumulates in float64. The published recipe's per-SNP sum is
# mathematically zero, so what it reports is rounding residue: float32 gives ~1e-5 that
# standardises into a convincing-looking scan, float64 gives ~1e-14 and an *uncorrelated*
# set of peaks. Its accumulator is therefore part of its answer, and it keeps float32 --
# the width every earlier run of this tool used -- so that "published" reproduces those
# outputs exactly. Nothing makes it agree with isoRelate or ibdutils, which are float64.
_DTYPE = np.float64
_VARIANT_DTYPE = {"corrected": np.float64, "published": np.float32}


def _variant_dtype(variant):
    return _DTYPE if variant == "corrected" else _VARIANT_DTYPE[variant]


def selection_bin_key(af, n_alleles=None, maf_k=None):
    """The frequency-bin key: ``maf`` at a biallelic site, ``1 - max(p_i)`` beyond one.

    Binning exists to compare a SNP against others of similar informativeness. Both keys
    order sites the same way in principle, so the obvious move is to bin on ``he`` -- and
    that is wrong for a reason only real data shows.

    **Heterozygosity is flat near p = 0.5.** MAFs of 0.48 and 0.52 have the *same* ``he``, so
    binning on it merges frequency classes that ``maf`` keeps apart. Synthetic frequencies
    drawn from a uniform are all distinct and hide this; a real panel is ``ac/an`` with ``an``
    around a thousand, so thousands of SNPs share each value. On the 27k-SNP Uganda callset,
    binning on ``he`` moved 248 SNPs between bins -- and since a bin's mean and sd shift with
    its membership, **12,986 z-scores changed**, by up to 4.6.

    So the key has to reduce to ``maf`` *exactly* at a biallelic site. ``1 - max(p_i)`` does:
    with two alleles it is the frequency of the less common one, which is what ``maf`` is. At
    k alleles it reaches ``1 - 1/k``, so a site more balanced than any biallelic one bins
    above them all, which is the improvement that was wanted.

    ``maf_k`` is that column from ``compute_allele_freqs``; without it, or at a site with two
    alleles, the ``maf`` computed here from ``af`` is used verbatim so the bins do not move.
    """
    maf = np.where(af <= 0.5, af, 1 - af)
    if maf_k is None or n_alleles is None:
        return maf
    use_k = (n_alleles > 2) & np.isfinite(maf_k)
    return np.where(use_k, maf_k, maf)


def _scale_and_gate(af, he):
    """The scale divisor and the validity gate, from heterozygosity when it is available.

    Allele frequency enters this statistic in exactly two places -- this divisor and the
    binning variable -- and both are heterozygosity terms wearing a biallelic disguise.
    ``p(1-p)`` is ``H/2`` where ``H = 1 - sum(p_i^2)``, so ``sqrt(he/2)`` is not an
    alternative to ``sqrt(af(1-af))`` but the same quantity written for any number of
    alleles. At a biallelic site the two are equal to the last bit, which is why supplying
    ``he`` moves no existing result.

    The gate changes with it, and this is the part that was quietly deleting data. ``af`` is
    the pooled non-reference frequency, so a site where the reference is absent altogether --
    4 C and 4 G, say -- collapses to ``af = 1.0`` and failed ``0 < af < 1``, though it is
    perfectly polymorphic and more informative than most of the panel. ``he > 0`` is the
    honest test for "nothing to see here", and a genuinely monomorphic site still fails it.
    """
    if he is None:
        valid = ~np.isnan(af) & (af > 0) & (af < 1)
        return valid, np.where(valid, np.sqrt(af * (1 - af)), np.nan)
    valid = ~np.isnan(he) & (he > 0)
    return valid, np.where(valid, np.sqrt(he / 2.0), np.nan)


def compute_selection_statistic(mat, af: np.ndarray, n_bins: int = 100,
                                label: str = "", variant: str = "corrected",
                                tail: str = "upper",
                                he: np.ndarray | None = None,
                                n_alleles: np.ndarray | None = None,
                                maf_k: np.ndarray | None = None) -> tuple[dict, pd.DataFrame]:
    """Per-SNP selection statistic. See the module docstring for ``variant`` and ``tail``.

    ``he`` is expected heterozygosity, ``1 - sum(p_i^2)`` over every allele including the
    reference (``compute_allele_freqs`` emits it). Supplying it generalises the two places
    allele frequency is used -- the scale divisor and the frequency binning -- to any number
    of alleles, and changes nothing on a biallelic panel. Without it the biallelic forms are
    used, which is what tables written before the column existed require.

    ``n_alleles`` and ``maf_k`` are the binning half of the same generalisation; see
    :func:`selection_bin_key` for why they are separate from ``he`` rather than derived from
    it. Without them the binning is unchanged.
    """
    if variant not in VARIANTS:
        raise ValueError(f"variant must be one of {VARIANTS}, got {variant!r}")
    n_pairs, n_snps = mat.shape
    tag = f"[{label}] " if label else ""
    print(f"  {tag}n_pairs={n_pairs:,}  n_snps={n_snps:,}")
    estimated_gb = n_pairs * n_snps * np.dtype(_variant_dtype(variant)).itemsize / 1e9
    if estimated_gb > 16:
        print(f"  {tag}~{estimated_gb:.1f} GB dense — using chunked path")
        return _compute_chunked(mat, af, n_bins, variant=variant, tail=tail, he=he,
                                n_alleles=n_alleles, maf_k=maf_k)
    print(f"  {tag}~{estimated_gb:.1f} GB dense — using dense path")
    return _compute_dense(mat, af, n_bins, variant=variant, tail=tail, he=he,
                          n_alleles=n_alleles, maf_k=maf_k)


def _compute_dense(mat, af, n_bins, variant="corrected", tail="upper", he=None,
                   n_alleles=None, maf_k=None):
    n_pairs, n_snps = mat.shape
    dt = _variant_dtype(variant)
    X = mat.toarray().astype(dt).T          # (snps, pairs)
    X -= X.mean(axis=0, keepdims=True)      # per-pair: removes each pair's relatedness
    if variant == "published":
        # Centring each SNP and then summing that same SNP cancels exactly; kept only to
        # reproduce isoRelate / ibdutils output.
        X -= X.mean(axis=1, keepdims=True)
    valid, denom = _scale_and_gate(af, he)
    denom = denom.astype(dt)
    X /= denom[:, np.newaxis]
    raw_stat = np.nansum(X, axis=1) / np.sqrt(n_pairs)
    raw_stat = np.where(valid, raw_stat, np.nan)
    return _normalise_and_finalise(raw_stat, af, valid, n_bins, tail=tail,
                                   n_alleles=n_alleles, maf_k=maf_k)


def _compute_chunked(mat, af, n_bins, chunk_size=500, variant="corrected", tail="upper",
                     he=None, n_alleles=None, maf_k=None):
    n_pairs, n_snps = mat.shape
    dt = _variant_dtype(variant)
    pair_means = np.asarray(mat.mean(axis=1)).ravel().astype(dt)
    valid, denom_all = _scale_and_gate(af, he)
    raw_stat = np.full(n_snps, np.nan, dtype=np.float64)
    n_chunks = (n_snps + chunk_size - 1) // chunk_size
    for c in range(n_chunks):
        start = c * chunk_size
        end = min(start + chunk_size, n_snps)
        if c % max(1, n_chunks // 10) == 0:
            print(f"    chunk {c+1}/{n_chunks}  SNPs {start}-{end}", end="\r", flush=True)
        chunk = mat[:, start:end].toarray().astype(dt).T  # (chunk, pairs)
        chunk -= pair_means[np.newaxis, :]
        if variant == "published":
            chunk -= chunk.mean(axis=1, keepdims=True)
        chunk_valid = valid[start:end]
        denom = denom_all[start:end].astype(dt)
        chunk /= denom[:, np.newaxis]
        rs = np.nansum(chunk, axis=1) / np.sqrt(n_pairs)
        raw_stat[start:end] = np.where(chunk_valid, rs, np.nan)
    print()
    return _normalise_and_finalise(raw_stat, af, valid, n_bins, tail=tail,
                                   n_alleles=n_alleles, maf_k=maf_k)


def _normalise_and_finalise(raw_stat, af, valid, n_bins, tail="upper",
                            n_alleles=None, maf_k=None):
    n_snps = len(raw_stat)
    maf = np.where(af <= 0.5, af, 1 - af)      # still reported, whatever the binning uses
    key = selection_bin_key(af, n_alleles=n_alleles, maf_k=maf_k)
    key_name = "maf" if maf_k is None or n_alleles is None else "binkey"

    bin_ids = np.full(n_snps, -1, dtype=int)
    valid_idx = np.where(valid)[0]
    key_valid = key[valid_idx]
    bin_edges = np.quantile(key_valid, np.linspace(0, 1, n_bins + 1))
    bin_edges[-1] += 1e-9
    bin_ids[valid_idx] = np.digitize(key_valid, bin_edges) - 1

    z_score = np.full(n_snps, np.nan)
    bin_records = []
    for b in range(n_bins):
        idx = np.where(bin_ids == b)[0]
        if len(idx) < 2:
            continue
        vals = raw_stat[idx]
        mu = np.nanmean(vals)
        sd = np.nanstd(vals, ddof=1)
        if sd == 0 or np.isnan(sd):
            continue
        z_score[idx] = (vals - mu) / sd
        bin_records.append({
            "bin": b, "n_snps": len(idx),
            f"{key_name}_min": bin_edges[b], f"{key_name}_max": bin_edges[b + 1],
            "mean": mu, "sd": sd,
        })
    bin_df = pd.DataFrame(bin_records)

    if tail not in TAILS:
        raise ValueError(f"tail must be one of {TAILS}, got {tail!r}")
    chi2_stat = z_score ** 2
    if tail == "two-sided":
        # chi2(1) on z^2 is the two-sided normal test, so a SNP shared LESS than expected
        # scores exactly like one shared more
        pval = np.where(~np.isnan(chi2_stat), chi2.sf(chi2_stat, df=1), np.nan)
    else:
        # excess sharing only: a deficit gets p -> 1 rather than a mirrored small p
        pval = np.where(~np.isnan(z_score), norm.sf(z_score), np.nan)
    neg_log10_p = np.where(pval > 0, -np.log10(pval), np.nan)
    direction = np.where(np.isnan(z_score), "",
                         np.where(z_score >= 0, "excess", "deficit"))

    return {
        "raw_stat": raw_stat, "z_score": z_score,
        "chi2_stat": chi2_stat, "pval": pval, "direction": direction,
        "neg_log10_p": neg_log10_p, "bin_id": bin_ids, "maf": maf,
    }, bin_df


def benjamini_hochberg(pval):
    """Benjamini-Hochberg q-values; `NaN` in, `NaN` out.

    BH controls the expected *proportion* of false positives among the SNPs called, where
    Bonferroni controls the probability of even one. It is valid under independence and
    under positive regression dependency -- the usual justification for using it on a
    genome scan, since SNPs in linkage are positively correlated. Treat that as an
    approximation: BH is now known not to control FDR in general for correlated two-sided
    tests. The larger practical worry is upstream of either correction, in whether the
    p-values are calibrated at all -- see :func:`genomic_inflation`.
    """
    p = np.asarray(pval, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.flatnonzero(~np.isnan(p))
    if not ok.size:
        return q
    m = ok.size
    order = ok[np.argsort(p[ok])]
    ranked = p[order] * m / np.arange(1, m + 1)
    q[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    return q


def genomic_inflation(chi2_stat):
    """Median chi2 over its null expectation -- 1.0 when the reference is right.

    Far from 1 means the chi2(1) null is the wrong distribution, and every p-value drawn
    from it is wrong with it. That happens easily here: the z-scores are standardised to
    unit variance within each MAF bin, which pins the first two moments and says nothing
    about the shape, while IBD sharing is autocorrelated and heavy-tailed. Read this
    before reading any threshold.
    """
    v = np.asarray(chi2_stat, dtype=float)
    v = v[~np.isnan(v)]
    if not v.size:
        return np.nan
    return float(np.median(v) / chi2.ppf(0.5, 1))


def permutation_null(mat, af, n_perm=200, n_bins=100, alpha=0.05, seed=0,
                     progress=None, variant="corrected", tail="upper"):
    """Null distribution for the selection statistic, built by moving the IBD around.

    Bonferroni and Benjamini-Hochberg both read their p-values off a chi2(1) that does not
    fit (see :func:`genomic_inflation`), and both treat SNPs as independent tests when one
    IBD segment spans hundreds of them. This generates the null rather than assuming it.

    Each replicate slides **every pair's IBD segments to a random position** along the SNP
    axis, wrapping at the end. That keeps each pair's total sharing, its segment count and
    its segment lengths exactly as observed -- so relatedness, block structure and the
    resulting autocorrelation all survive -- and destroys only the thing being tested:
    whether pairs share *the same* locus. The statistic is then recomputed from scratch,
    MAF binning and all. Expect the result to be stricter than the parametric lines, and
    possibly far stricter: where the replicates' genome-wide maxima sit above the Bonferroni
    line, that line's true family-wise error rate is not the nominal 5%. Compare the two
    rather than assuming a margin.

    Four summaries come out of the one pass, trading resolution against assumptions:

    ``threshold``
        The ``1 - alpha`` quantile of the per-replicate genome-wide **maxima**: the largest
        score reachable with no locus-specific sharing. Family-wise, and the one summary
        that pools nothing across SNPs.
    ``p_pointwise``
        Per SNP, how often the null at *that same SNP* reached the observed value. Assumes
        nothing beyond the shift, but cannot resolve below ``1 / (n_perm + 1)``, so it is a
        check on the top hits rather than an input to FDR.
    ``p_stratified``
        The same count within the SNP's own MAF bin. Also assumption-free, and ``n_bins``
        times finer than pointwise -- but the floor is per-SNP, and MAF ties make the bins
        very unequal, so SNPs in small bins can be unreachable at any score.
    ``p_pooled``
        Per SNP, how often *any* null value anywhere reached the observed value. Resolution
        ``1 / (n_perm * n_snps + 1)``, fine enough to feed Benjamini-Hochberg, at the price
        of assuming the null is exchangeable across MAF bins. ``bin_tail_rate`` reports
        whether it is.

    All three p-values use the Phipson-Smyth ``(1 + exceedances) / (1 + draws)`` form, so
    none is ever 0 -- a permutation cannot evidence a p below its own resolution.

    ``bin_tail_rate`` is the exchangeability check behind ``p_pooled``: the share of each
    MAF bin's nulls above one common reference (the 99th percentile of the first
    replicate). Exchangeability puts every bin at 0.01; a bin far off has a differently
    shaped null, and pooling across it makes its SNPs' p-values too small. Bins with fewer
    than 1000 draws are left ``NaN``, since a 1% rate is unresolvable below that.

    Args:
        mat: The (pairs x SNPs) binary matrix for one group.
        af: Allele frequencies aligned to its columns.
        n_perm: Replicates. 200 suffices for a 5% quantile. Benjamini-Hochberg over
            ``p_pooled`` needs ``n_perm >= 10 / q`` for an order of magnitude of headroom,
            because the smallest reachable q is about ``1 / n_perm`` whatever the cohort
            size.
        n_bins, alpha, seed: MAF bins, family-wise level, RNG seed.
        progress: Optional ``callable(i, n_perm)`` for a progress line.

    Returns:
        A dict with ``threshold``, ``maxima``, ``p_pointwise``, ``p_pooled``,
        ``p_stratified``, ``n_stratified``, ``bin_tail_rate``, ``bin_tail_rate_skipped``,
        ``n_perm`` and ``n_pool``. Everything per-SNP is aligned to the columns of ``mat``
        and carries ``NaN`` wherever the observed statistic does.
    """
    from scipy import sparse

    coo = mat.tocoo()
    n_pairs, n_snps = mat.shape
    obs = compute_selection_statistic(mat, af, n_bins=n_bins, label="", variant=variant,
                                      tail=tail)[0]["neg_log10_p"]
    obs_ok = np.isfinite(obs)

    rng = np.random.default_rng(seed)
    maxima = np.empty(n_perm)
    point = np.zeros(n_snps, dtype=np.int64)     # null beat the observed AT this SNP
    pooled = np.zeros(n_snps, dtype=np.int64)    # null beat it ANYWHERE
    strat = np.zeros(n_snps, dtype=np.int64)     # ...anywhere IN ITS OWN MAF BIN
    bin_n = np.zeros(n_bins, dtype=np.int64)
    bin_hi = np.zeros(n_bins, dtype=np.int64)
    n_pool = 0
    ref = np.nan
    # Bin membership is a function of `af` and `n_bins` alone, so it is the same in every
    # replicate as in the observed scan -- which is what lets the stratified count below
    # compare a SNP only against nulls from its own bin.
    members = None

    for r in range(n_perm):
        shift = rng.integers(0, n_snps, size=n_pairs)
        col = (coo.col + shift[coo.row]) % n_snps
        null = sparse.coo_matrix((coo.data, (coo.row, col)), shape=mat.shape).tocsr()
        st, _ = compute_selection_statistic(null, af, n_bins=n_bins, label="",
                                            variant=variant, tail=tail)
        v = st["neg_log10_p"]
        maxima[r] = np.nanmax(v)

        # NaN compares False, so unusable SNPs simply never count as exceedances
        point += v >= obs
        good = np.sort(v[np.isfinite(v)])
        n_pool += good.size
        # searchsorted puts NaN past the end, giving those SNPs a count of 0; they are
        # masked to NaN below rather than silently reported as significant
        pooled += good.size - np.searchsorted(good, obs, side="left")

        b = st["bin_id"]
        keep = np.isfinite(v) & (b >= 0)
        if r == 0:
            ref = float(np.nanquantile(v, 0.99))
            members = [np.where(b == k)[0] for k in range(n_bins)]
        bin_n += np.bincount(b[keep], minlength=n_bins)[:n_bins]
        bin_hi += np.bincount(b[keep & (v >= ref)], minlength=n_bins)[:n_bins]

        for idx in members:
            if idx.size == 0:
                continue
            vb = np.sort(v[idx][np.isfinite(v[idx])])
            if vb.size:
                strat[idx] += vb.size - np.searchsorted(vb, obs[idx], side="left")

        if progress is not None:
            progress(r + 1, n_perm)

    # draws behind each SNP's stratified p: its own bin's finite nulls, over all replicates
    n_strat = np.zeros(n_snps, dtype=np.int64)
    for k, idx in enumerate(members or []):
        n_strat[idx] = bin_n[k]

    p_point = np.where(obs_ok, (1 + point) / (1 + n_perm), np.nan)
    p_pool = np.where(obs_ok, (1 + pooled) / (1 + n_pool), np.nan)
    p_strat = np.where(obs_ok & (n_strat > 0),
                       (1 + strat) / (1 + np.maximum(n_strat, 1)), np.nan)
    # Only bins with enough draws to resolve a 1% rate get one: below ~1000 the estimate
    # is dominated by counting noise, and a small bin returning 0 would otherwise make the
    # spread look infinite when nothing is actually wrong with it.
    _MIN_DRAWS = 1000
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = np.where(bin_n >= _MIN_DRAWS, bin_hi / np.maximum(bin_n, 1), np.nan)

    return {
        "threshold": float(np.quantile(maxima, 1 - alpha)),
        "maxima": maxima,
        "p_pointwise": p_point,
        "p_pooled": p_pool,
        "p_stratified": p_strat,
        "n_stratified": n_strat,
        "bin_tail_rate": rate,
        "bin_tail_rate_skipped": int(np.sum((bin_n > 0) & (bin_n < _MIN_DRAWS))),
        "n_perm": int(n_perm),
        "n_pool": int(n_pool),
    }


def assemble_output(snp_df, stats, alpha, group=None, fdr_alpha=None, perm=None,
                    pool="global", variant="corrected", tail="upper"):
    """Per-SNP table plus every threshold and the calibration diagnostic.

    ``perm`` is the dict from :func:`permutation_null`, or ``None`` to skip the
    permutation columns entirely. ``pool`` chooses which of its empirical p-values drives
    ``q_empirical``: ``"global"`` pools the null across all SNPs (fine resolution, assumes
    exchangeability) and ``"bin"`` stays inside each MAF bin (assumption-free, ``n_bins``
    times coarser). Both are written either way; only the q-value follows ``pool``.

    Returns ``(df, info)``. ``info`` holds each cutoff on the ``-log10(p)`` scale so any
    of them can be drawn as a line, and ``lambda_gc``, which says how much the parametric
    ones are worth.
    """
    fdr_alpha = alpha if fdr_alpha is None else fdr_alpha
    n_valid = int(np.sum(~np.isnan(stats["chi2_stat"])))
    threshold = -np.log10(alpha / n_valid) if n_valid > 0 else np.nan

    out = snp_df.copy()
    if group is not None:
        out.insert(0, "group", group)
    out["maf"] = stats["maf"]
    out["bin_id"] = stats["bin_id"]
    out["raw_stat"] = stats["raw_stat"]
    out["z_score"] = stats["z_score"]
    out["chi2_stat"] = stats["chi2_stat"]
    out["pval"] = stats["pval"]
    if "direction" in stats:
        out["direction"] = stats["direction"]
    out["neg_log10_p"] = stats["neg_log10_p"]
    out["q_value"] = benjamini_hochberg(stats["pval"])
    out["significant"] = out["neg_log10_p"] >= threshold          # Bonferroni (FWER)
    out["significant_fdr"] = out["q_value"] < fdr_alpha           # Benjamini-Hochberg
    perm_threshold = perm["threshold"] if perm else None
    if perm_threshold is not None and np.isfinite(perm_threshold):
        out["significant_perm"] = out["neg_log10_p"] >= perm_threshold
    if perm:
        # Calibrated companions to `pval`/`q_value`: same SNPs, same order, but a
        # reference drawn from the data rather than assumed. `p_pointwise` bottoms out at
        # 1/(n_perm+1) so it cannot drive FDR; it is here as the assumption-free check on
        # whatever the pooled p-value calls.
        if pool not in ("global", "bin"):
            raise ValueError("pool must be 'global' or 'bin'")
        out["p_pointwise"] = perm["p_pointwise"]
        out["p_empirical"] = perm["p_pooled"]
        out["p_empirical_binned"] = perm["p_stratified"]
        chosen = perm["p_pooled"] if pool == "global" else perm["p_stratified"]
        out["q_empirical"] = benjamini_hochberg(chosen)
        out["significant_fdr_perm"] = out["q_empirical"] < fdr_alpha

    # The BH critical value for the number of rejections, k * q / m, as -log10(p) so it
    # plots as a line. Not the largest p actually called: this is the cutoff BH applies,
    # so `neg_log10_p >= line` reproduces the flag exactly. It is always at or below the
    # Bonferroni line (which is the k = 1 case).
    n_rej = int(out["significant_fdr"].fillna(False).sum())
    fdr_threshold = (float(-np.log10(n_rej * fdr_alpha / n_valid))
                     if n_rej and n_valid else np.nan)

    # The smallest q the permutation could possibly produce: the top SNP's p bottoms out
    # at 1/(1 + n_perm * n_snps), and BH multiplies the rank-1 p-value by m. So the whole
    # empirical-FDR column is dead unless this sits below `fdr_alpha` -- and it wants to
    # sit well below, since one stray null exceedance at the top SNP multiplies it. In
    # round terms q_floor ~ 1 / n_perm, so n_perm >= 10 / fdr_alpha buys an order of
    # magnitude of headroom (200 replicates at the default q < 0.05).
    # The smallest q a *lone* extreme SNP could reach: BH gives rank k the value p*m/k, and
    # p bottoms out at 1/(1 + draws behind it), so at rank 1 that is m/(1 + draws). A block
    # of k SNPs tied at the resolution limit reaches k times lower, which is why this is a
    # conservative bound rather than a hard cutoff. Under global pooling every SNP has the
    # same draws and so the same bound; under bin pooling it is per-SNP, and MAF bins are
    # wildly unequal in practice (ties collapse the quantile edges), so a single "best case"
    # number would hide most of the column being out of reach.
    # How the MAF binning actually came out. Requesting N equal-frequency bins does not
    # give N: allele frequency is k/n for a smallish n, so MAF ties collapse the quantile
    # edges. This is a property of the statistic, not of the permutation, so it is
    # computed either way -- the within-bin standardisation is step 5 of the method.
    usable = ~np.isnan(stats["chi2_stat"])
    counts = np.bincount(stats["bin_id"][usable & (stats["bin_id"] >= 0)])
    n_bins_used = int((counts > 0).sum())
    big_bin = float(counts.max() / n_valid) if counts.size and n_valid else np.nan

    q_floor, frac_dead = np.nan, np.nan
    if perm and n_valid:
        draws = (np.full(len(stats["bin_id"]), perm["n_pool"]) if pool == "global"
                 else perm["n_stratified"])
        with np.errstate(invalid="ignore", divide="ignore"):
            per_snp = np.where(draws > 0, n_valid / (1.0 + draws), np.inf)
        q_floor = float(np.min(per_snp[usable]))
        frac_dead = float(np.mean(per_snp[usable] > fdr_alpha))

    # The empirical-FDR line, as the smallest observed score BH kept. Unlike the k*q/m
    # form above there is no closed expression for it, because the empirical p-values are
    # a step function of the score rather than a smooth transform of it.
    emp_line = np.nan
    if perm and out["significant_fdr_perm"].any():
        emp_line = float(out.loc[out["significant_fdr_perm"].fillna(False),
                                 "neg_log10_p"].min())

    info = {
        "alpha": alpha, "n_tests": n_valid, "neg_log10_p_threshold": threshold,
        "neg_log10_p_perm_threshold": (float(perm_threshold)
                                       if perm_threshold is not None else np.nan),
        "n_significant_perm": (int(out["significant_perm"].sum())
                               if "significant_perm" in out else -1),
        "fdr_alpha": fdr_alpha, "neg_log10_p_fdr_threshold": fdr_threshold,
        "n_significant": int(out["significant"].sum()),
        "n_significant_fdr": int(out["significant_fdr"].fillna(False).sum()),
        "neg_log10_p_emp_fdr_threshold": emp_line,
        "n_significant_fdr_perm": (int(out["significant_fdr_perm"].fillna(False).sum())
                                   if perm else -1),
        "n_perm": perm["n_perm"] if perm else 0,
        "p_empirical_resolution": (1.0 / (1 + perm["n_pool"])) if perm else np.nan,
        "q_empirical_floor": q_floor, "empirical_pool": (pool if perm else ""),
        "xirs_variant": variant, "tail": tail,
        "frac_q_unreachable": frac_dead,
        "n_bins_used": n_bins_used, "largest_bin_frac": big_bin,
        "perm_bin_tail_min": (float(np.nanmin(perm["bin_tail_rate"])) if perm else np.nan),
        "perm_bin_tail_max": (float(np.nanmax(perm["bin_tail_rate"])) if perm else np.nan),
        "lambda_gc": genomic_inflation(stats["chi2_stat"]),
    }
    return out, info
