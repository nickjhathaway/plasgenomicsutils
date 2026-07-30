"""Within-sample AD cleaning + re-genotyping over a whole VCF/BCF.

A site's ``DP`` can look deep while the depth actually used to genotype (``AD`` /
its sum ``ADS``) is tiny — such sites are almost always repetitive/artifact
regions. This judges depth by AD, zeros sub-threshold allele depths per sample,
recomputes ADS, and re-genotypes conservatively from the cleaned counts so a
little genotyping error does not force discarding otherwise-good biallelic SNPs.

Uses cyvcf2 for bulk per-record FORMAT-array access (read AD/ADS and write
AD/ADS/GT as whole (n_samples x n_alleles) arrays), which is far faster than
per-sample access on large cohorts. The re-genotyping matches
:func:`plasgenomicsutils.lib.ad_genotype.regenotype_from_ad` exactly.
"""

from __future__ import annotations

import numpy as np

from .ad_genotype import clean_ad_matrix, regenotype_from_ad, regenotype_matrix

_MISS = np.iinfo(np.int32).min  # cyvcf2 sentinel for a missing integer FORMAT value


def filter_ad_regenotype(input_vcf: str, output_vcf: str, *, min_reads: int = 2,
                         min_freq: float = 0.01, het_min_af: float = 0.2,
                         restrict_to_called: bool = False) -> None:
    """Clean AD and re-genotype every record; records lacking AD/ADS pass through.

    The frequency denominator is the ``ADS`` FORMAT field (summed genotyping
    depth). ADS is recomputed from the cleaned AD before writing. A record is
    passed through unchanged if it lacks AD or ADS, or any sample's AD/ADS is
    missing.

    ``restrict_to_called`` limits each sample's genotype to the alleles the
    upstream caller already used, so calls are only narrowed, never promoted to a
    new allele. Use it when the caller's likelihood-based genotypes should be
    trusted over raw read counts.
    """
    from cyvcf2 import VCF, Writer

    vcf = VCF(input_vcf)
    out = Writer(output_vcf, vcf)
    for v in vcf:
        ad = v.format("AD")
        ads = v.format("ADS")
        if ad is None or ads is None:
            out.write_record(v)
            continue
        ad = ad.astype(np.int64)
        depth = ads[:, 0].astype(np.int64)
        if (ad < 0).any() or (depth < 0).any():  # any missing AD/ADS -> pass through
            out.write_record(v)
            continue

        cleaned = clean_ad_matrix(ad.astype(float), depth.astype(float),
                                  min_reads, min_freq, protect_ref=False)
        new_ads = cleaned.sum(axis=1).astype(np.int32)
        gts = v.genotypes  # list of [a1, a2, ..., phased]

        if restrict_to_called:
            for i in range(len(gts)):
                called = tuple(gts[i][:-1])
                g = regenotype_from_ad([int(x) for x in cleaned[i]], het_min_af,
                                       called_alleles=called)
                gts[i] = [-1, -1, gts[i][-1]] if g is None else [g[0], g[1], gts[i][-1]]
        else:
            gt_a, gt_b = regenotype_matrix(cleaned, het_min_af)
            for i in range(len(gts)):
                gts[i] = [int(gt_a[i]), int(gt_b[i]), gts[i][-1]]

        v.set_format("AD", cleaned.astype(np.int32))
        v.set_format("ADS", new_ads)
        v.genotypes = gts
        out.write_record(v)

    out.close()
    vcf.close()
