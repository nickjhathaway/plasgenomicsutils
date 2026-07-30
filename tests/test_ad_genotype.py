"""Tests for the shared AD cleaning + re-genotyping core."""

import numpy as np

from plasgenomicsutils.lib.ad_genotype import clean_ad_matrix, regenotype_from_ad


def test_regenotype_missing_and_hom():
    assert regenotype_from_ad([0, 0]) is None
    assert regenotype_from_ad([10, 0]) == (0, 0)
    assert regenotype_from_ad([0, 10]) == (1, 1)


def test_regenotype_het_threshold():
    # minor AF 0.1 < 0.2 -> homozygous major (ref)
    assert regenotype_from_ad([9, 1]) == (0, 0)
    # minor AF 0.5 -> heterozygous
    assert regenotype_from_ad([5, 5]) == (0, 1)
    # minor AF ~0.23 >= 0.2 -> heterozygous
    assert regenotype_from_ad([10, 3]) == (0, 1)


def test_regenotype_multiallelic_and_custom_threshold():
    # two alts tie -> het of the top two allele indices (sorted)
    assert regenotype_from_ad([0, 5, 5]) == (1, 2)
    # a stricter threshold turns a 0.25 minor into homozygous major
    assert regenotype_from_ad([9, 3], het_min_af=0.3) == (0, 0)


def test_clean_ad_matrix_zeros_subthreshold():
    ad = np.array([[10, 1], [5, 5]], dtype=float)
    depth = np.array([11, 10], dtype=float)
    out = clean_ad_matrix(ad, depth, min_reads=2, min_freq=0.01, protect_ref=False)
    # row0 alt has 1 read (< 2) -> zeroed; row1 both kept
    assert out.tolist() == [[10, 0], [5, 5]]


def test_clean_ad_matrix_protect_ref():
    ad = np.array([[1, 20]], dtype=float)   # ref has only 1 read
    depth = np.array([21], dtype=float)
    protected = clean_ad_matrix(ad, depth, min_reads=2, min_freq=0.01, protect_ref=True)
    assert protected.tolist() == [[1, 20]]   # ref never zeroed
    unprotected = clean_ad_matrix(ad, depth, min_reads=2, min_freq=0.01, protect_ref=False)
    assert unprotected.tolist() == [[0, 20]]  # ref zeroed (1 < 2)
