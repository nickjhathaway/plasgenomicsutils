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

import re
import sys
from math import comb

import numpy as np

from .ad_genotype import clean_ad_matrix, regenotype_from_ad, regenotype_matrix

_MISS = np.iinfo(np.int32).min  # cyvcf2 sentinel for a missing integer FORMAT value

# genotype-linked (Number=G) FORMAT fields whose likelihoods are meaningless once we
# re-genotype from AD. Re-genotyping writes a diploid GT, so leaving a caller's PL from a
# different ploidy (e.g. hexaploid calls forced to diploid) makes the record internally
# inconsistent and breaks `bcftools view --trim-alt-alleles` downstream. We blank them to
# a diploid-consistent length on every re-genotyped record.
_GT_LIKELIHOOD_TAGS = ("PL", "GL", "GP", "PP", "PG")


def _present_likelihood_tags(vcf) -> list[str]:
    present = []
    for t in _GT_LIKELIHOOD_TAGS:
        try:
            if vcf.get_header_type(t) is not None:
                present.append(t)
        except KeyError:
            pass
    return present


def _format_ids(vcf) -> set:
    """The FORMAT tag IDs a cyvcf2 VCF's header declares."""
    return set(re.findall(r"^##FORMAT=<ID=([^,>]+)", vcf.raw_header, flags=re.M))


def filter_ad_regenotype(input_vcf: str, output_vcf: str, *, min_reads: int = 2,
                         min_freq: float = 0.01, het_min_af: float = 0.2,
                         restrict_to_called: bool = False,
                         drop_stale_likelihoods: bool = True,
                         ploidy: int | None = None,
                         add_ads: bool = True) -> None:
    """Clean AD and re-genotype every record; records lacking AD pass through.

    The frequency denominator is the ``ADS`` FORMAT field (summed genotyping
    depth). ADS is recomputed from the cleaned AD before writing. A record is
    passed through unchanged if it lacks AD, or any sample's AD/ADS is missing.

    ``add_ads`` (default ``True``) derives ADS from AD -- their sum, which is what
    ``singleton_filter_add_ads`` computes -- when the callset does not carry it, and
    adds the header line. Without this a callset with no ADS is passed through
    untouched: every record silently unfiltered, which looks exactly like a filter
    that had nothing to do. Set it ``False`` to keep that pass-through behaviour.

    ``restrict_to_called`` limits each sample's genotype to the alleles the
    upstream caller already used, so calls are only narrowed, never promoted to a
    new allele. Use it when the caller's likelihood-based genotypes should be
    trusted over raw read counts.

    ``drop_stale_likelihoods`` (default ``True``) blanks genotype-linked FORMAT
    fields (``PL``/``GL``/…) on every re-genotyped record: writing the output GT
    makes a caller's original-ploidy likelihoods inconsistent, which would break
    ``bcftools view --trim-alt-alleles`` later. They are unused downstream.

    ``ploidy`` sets the output genotype ploidy. ``None`` (default) keeps the
    conventional diploid coding used for *Plasmodium* (``0/1`` = mixed infection).
    When given it is validated against the input ploidy per record: a request
    **greater** than the input ploidy raises an error (genotypes/likelihoods cannot
    be promoted); **equal** is fine; **less** is allowed but warns and trims the
    genotype-linked fields to the lower ploidy. ``1`` calls the single
    best-supported allele (haploid); ``2`` is diploid.
    """
    from cyvcf2 import VCF, Writer

    out_ploidy = 2 if ploidy is None else int(ploidy)
    if out_ploidy not in (1, 2):
        raise SystemExit("ERROR: --ploidy must be 1 (haploid) or 2 (diploid)")

    vcf = VCF(input_vcf)
    n_samples = len(vcf.samples)
    gl_tags = _present_likelihood_tags(vcf) if drop_stale_likelihoods else []
    # asked before the header is touched: cyvcf2 raises KeyError on v.format() for a tag
    # the header does not declare, so this has to gate every read of it
    in_ads = "ADS" in _format_ids(vcf)
    derive_ads = add_ads and not in_ads
    if derive_ads:
        vcf.add_format_to_header({
            "ID": "ADS", "Number": "1", "Type": "Integer",
            "Description": "Summed allelic depth (sum of AD), the reads used to genotype"})
        sys.stderr.write("NOTE: no FORMAT/ADS in the input; deriving it as the sum of AD\n")
    out = Writer(output_vcf, vcf)
    warned_reduce = False
    for v in vcf:
        ad = v.format("AD")
        if ad is None:
            out.write_record(v)
            continue
        ad = ad.astype(np.int64)
        ads = v.format("ADS") if in_ads else None
        if ads is None:
            if not derive_ads:
                out.write_record(v)
                continue
            # the sum of this sample's own AD, which is what ADS means
            depth = np.where(ad < 0, 0, ad).sum(axis=1).astype(np.int64)
        else:
            depth = ads[:, 0].astype(np.int64)
        if (ad < 0).any() or (depth < 0).any():  # any missing AD/ADS -> pass through
            out.write_record(v)
            continue

        cleaned = clean_ad_matrix(ad.astype(float), depth.astype(float),
                                  min_reads, min_freq, protect_ref=False)
        new_ads = cleaned.sum(axis=1).astype(np.int32)
        gts = v.genotypes  # list of [a1, a2, ..., phased]

        if ploidy is not None:
            in_ploidy = max((len(g) - 1 for g in gts), default=out_ploidy)
            if out_ploidy > in_ploidy:
                raise SystemExit(
                    f"ERROR: --ploidy {out_ploidy} exceeds the input genotype ploidy "
                    f"({in_ploidy}) at {v.CHROM}:{v.POS}; genotypes/likelihoods cannot be "
                    "promoted to a higher ploidy.")
            if out_ploidy < in_ploidy and not warned_reduce:
                sys.stderr.write(
                    f"WARNING: reducing genotype ploidy {in_ploidy} -> {out_ploidy}; "
                    "genotype-linked FORMAT fields (PL/GL/…) are trimmed to match.\n")
                warned_reduce = True

        if out_ploidy == 1:
            major = np.argmax(cleaned, axis=1)          # ties -> lowest allele index
            totals = cleaned.sum(axis=1)
            for i in range(len(gts)):
                a = -1 if totals[i] == 0 else int(major[i])
                gts[i] = [a, gts[i][-1]]
        elif restrict_to_called:
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
        if gl_tags:
            # Number=G length for this record at the output ploidy
            n_alleles = len(v.ALT) + 1
            g_width = comb(n_alleles - 1 + out_ploidy, out_ploidy)
            miss = np.full((n_samples, g_width), _MISS, dtype=np.int32)
            for t in gl_tags:
                if v.format(t) is not None:
                    v.set_format(t, miss)
        out.write_record(v)

    out.close()
    vcf.close()
