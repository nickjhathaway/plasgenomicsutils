"""Within-sample AD cleaning + re-genotyping over a whole VCF/BCF (streaming).

A site's ``DP`` can look deep while the depth actually used to genotype (``AD`` /
its sum ``ADS``) is tiny — such sites are almost always repetitive/artifact
regions. This judges depth by AD, zeros sub-threshold allele depths per sample,
recomputes ADS, and re-genotypes conservatively from the cleaned counts so a
little genotyping error does not force discarding otherwise-good biallelic SNPs.
"""

from __future__ import annotations

import numpy as np

from .ad_genotype import clean_ad_matrix, regenotype_from_ad


def filter_ad_regenotype(input_vcf: str, output_vcf: str, *, min_reads: int = 2,
                         min_freq: float = 0.01, het_min_af: float = 0.2,
                         restrict_to_called: bool = False) -> None:
    """Clean AD and re-genotype every record; records lacking AD/ADS pass through.

    The frequency denominator is the ``ADS`` FORMAT field (summed genotyping
    depth). ADS is recomputed from the cleaned AD before writing.

    ``restrict_to_called`` limits each sample's genotype to the alleles the
    upstream caller already used, so calls are only narrowed, never promoted to a
    new allele. Use it when the caller's likelihood-based genotypes should be
    trusted over raw read counts.
    """
    import pysam

    vcf = pysam.VariantFile(input_vcf)
    out = pysam.VariantFile(output_vcf, _write_mode(output_vcf), header=vcf.header)

    snames = list(vcf.header.samples)
    for rec in vcf:
        n_alleles = len(rec.alleles)
        ad_rows, ads_vals, ok = _gather(rec, snames, n_alleles)
        if not ok:
            out.write(rec)
            continue

        cleaned = clean_ad_matrix(ad_rows, ads_vals, min_reads, min_freq, protect_ref=False)
        new_ads = cleaned.sum(axis=1).astype(int)

        for i, sname in enumerate(snames):
            sample = rec.samples[sname]
            row = [int(x) for x in cleaned[i]]
            called = tuple(sample["GT"]) if restrict_to_called else None
            sample["AD"] = tuple(row)
            sample["ADS"] = int(new_ads[i])
            phased = bool(sample.phased)
            gt = regenotype_from_ad(row, het_min_af, called_alleles=called)
            sample["GT"] = (None, None) if gt is None else gt
            sample.phased = phased
        out.write(rec)

    out.close()
    vcf.close()


def _gather(rec, snames, n_alleles):
    """Collect per-sample AD (n_samples x n_alleles) and ADS depth denominators."""
    ad_rows = np.zeros((len(snames), n_alleles), dtype=float)
    ads_vals = np.zeros(len(snames), dtype=float)
    for i, sname in enumerate(snames):
        sample = rec.samples[sname]
        ad = sample.get("AD", None)
        ads = sample.get("ADS", None)
        if ad is None or ads is None or any(v is None for v in ad) or len(ad) != n_alleles:
            return None, None, False
        ad_rows[i, :] = ad
        ads_vals[i] = ads[0] if isinstance(ads, (tuple, list)) else ads
    return ad_rows, ads_vals, True


def _write_mode(path: str) -> str:
    if path.endswith(".bcf"):
        return "wb"
    if path.endswith(".vcf.gz"):
        return "wz"
    return "w"
