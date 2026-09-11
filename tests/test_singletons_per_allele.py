"""Singletons are a property of an allele, not of a record.

"Exactly one sample carries the non-reference allele" and "exactly one sample carries *a*
non-reference allele" are the same question at a biallelic site and different questions as
soon as there are two alternates. At a triallelic site where A carries ALT1 and B carries
ALT2, the truth is two singletons; counting carriers per record saw two non-reference samples
and booked a **doubleton shared by A and B**.

That was wrong twice over for what the module is for: it deflated both samples' singleton
counts, and it fed `shared[(i, j)]`, which drives the "near-identical to X" flag. A cohort
with many multiallelic sites accumulated pairwise evidence of near-identity between samples
carrying *different* alleles.
"""

import pytest

pytest.importorskip("cyvcf2")

from plasgenomicsutils.lib.singletons import count_singletons

_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
)


def _vcf(tmp_path, samples, rows, name="in.vcf"):
    hdr = _HDR + ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(samples) + "\n")
    body = "".join(f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts) + "\n"
                   for pos, alt, gts in rows)
    p = tmp_path / name
    p.write_text(hdr + body)
    return str(p)


def test_two_alternates_private_to_different_samples_are_two_singletons(tmp_path):
    samples = ["s1", "s2", "s3", "s4"]
    rows = [(1000, "C,G", ["1/1", "2/2", "0/0", "0/0"])]
    df, n = count_singletons(_vcf(tmp_path, samples, rows), max_missing_frac=0.5)
    got = df.set_index("sample")["n_singleton"].to_dict()
    assert got["s1"] == 1 and got["s2"] == 1
    assert got["s3"] == 0 and got["s4"] == 0


def test_they_are_not_recorded_as_sharing_a_doubleton(tmp_path):
    """The `near-identical to X` flag is built on this, so a fabricated pair is expensive."""
    samples = ["s1", "s2", "s3", "s4"]
    rows = [(1000 + i, "C,G", ["1/1", "2/2", "0/0", "0/0"]) for i in range(20)]
    df, _n = count_singletons(_vcf(tmp_path, samples, rows), max_missing_frac=0.5)
    d = df.set_index("sample")
    assert d.loc["s1", "n_doubleton"] == 0
    assert d.loc["s2", "n_doubleton"] == 0


def test_one_alternate_carried_by_two_samples_is_still_a_doubleton(tmp_path):
    samples = ["s1", "s2", "s3", "s4"]
    rows = [(1000, "C,G", ["1/1", "1/1", "0/0", "0/0"])]
    df, _n = count_singletons(_vcf(tmp_path, samples, rows), max_missing_frac=0.5)
    d = df.set_index("sample")
    assert d.loc["s1", "n_doubleton"] == 1 and d.loc["s2", "n_doubleton"] == 1
    assert d.loc["s1", "n_singleton"] == 0
    assert d.loc["s1", "top_partner"] == "s2"


def test_a_het_carrying_both_alternates_is_a_carrier_of_each(tmp_path):
    """`1/2` is one sample carrying two alleles, so each allele has one carrier here."""
    samples = ["s1", "s2", "s3"]
    rows = [(1000, "C,G", ["1/2", "0/0", "0/0"])]
    df, _n = count_singletons(_vcf(tmp_path, samples, rows), max_missing_frac=0.5)
    d = df.set_index("sample")
    assert d.loc["s1", "n_singleton"] == 2      # private in both C and G
    assert d.loc["s2", "n_singleton"] == 0


def test_a_biallelic_cohort_counts_exactly_as_before(tmp_path):
    samples = ["s1", "s2", "s3", "s4"]
    rows = [
        (1000, "C", ["1/1", "0/0", "0/0", "0/0"]),   # s1 singleton
        (2000, "C", ["1/1", "1/1", "0/0", "0/0"]),   # s1+s2 doubleton
        (3000, "C", ["0/0", "0/0", "1/1", "0/0"]),   # s3 singleton
        (4000, "C", ["1/1", "1/1", "1/1", "0/0"]),   # neither
    ]
    df, n = count_singletons(_vcf(tmp_path, samples, rows), max_missing_frac=0.5)
    d = df.set_index("sample")
    assert n == 4
    assert list(d.loc[["s1", "s2", "s3", "s4"], "n_singleton"]) == [1, 0, 1, 0]
    assert list(d.loc[["s1", "s2", "s3", "s4"], "n_doubleton"]) == [1, 1, 0, 0]
    assert d.loc["s1", "top_partner"] == "s2"


# --- `*` is not an allele -------------------------------------------------------------
#
# A spanning deletion says the sequence is absent on that haplotype. Counting it as an
# alternate does two kinds of damage here. The first is a private "SNP" for a sample that
# carries no SNP. The second is worse: two samples that share a deletion are booked a
# DOUBLETON, and doubletons feed `shared[(i, j)]`, which drives the "near-identical to X"
# flag. Deletions are haplotype markers, so sharing one is common, and the fabricated
# evidence accumulates between samples that are not near-identical at all.


def _vcf_star(tmp_path, samples, rows, name="star.vcf"):
    hdr = _HDR + ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(samples) + "\n")
    body = "".join(f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts) + "\n"
                   for pos, alt, gts in rows)
    p = tmp_path / name
    p.write_text(hdr + body)
    return str(p)


def test_a_lone_spanning_deletion_is_not_a_singleton(tmp_path):
    samples = ["s1", "s2", "s3", "s4"]
    # s1 is the only sample with the deletion, and the only sample carrying the real ALT
    rows = [(1000, "T,*", ["2/2", "1/1", "0/0", "0/0"])]
    df, _n = count_singletons(_vcf_star(tmp_path, samples, rows), max_missing_frac=0.5)
    got = df.set_index("sample")["n_singleton"].to_dict()
    assert got["s1"] == 0        # the deletion is not a private variant
    assert got["s2"] == 1        # the real ALT still is


def test_two_samples_sharing_a_deletion_are_not_near_identical(tmp_path):
    samples = ["s1", "s2", "s3", "s4"]
    rows = [(1000 + i, "T,*", ["2/2", "2/2", "0/0", "0/0"]) for i in range(20)]
    df, _n = count_singletons(_vcf_star(tmp_path, samples, rows), max_missing_frac=0.5)
    d = df.set_index("sample")
    assert d.loc["s1", "n_doubleton"] == 0
    assert d.loc["s2", "n_doubleton"] == 0


def test_a_record_whose_only_alternate_is_a_deletion_is_skipped(tmp_path):
    samples = ["s1", "s2", "s3", "s4"]
    rows = [(1000, "*", ["1/1", "0/0", "0/0", "0/0"])]
    df, _n = count_singletons(_vcf_star(tmp_path, samples, rows), max_missing_frac=0.5)
    assert df["n_singleton"].sum() == 0
    assert df.attrs["n_star_only_records"] == 1


def test_the_real_alternates_beside_a_deletion_are_still_counted(tmp_path):
    samples = ["s1", "s2", "s3", "s4"]
    rows = [(1000, "T,*,G", ["1/1", "2/2", "3/3", "0/0"])]
    df, _n = count_singletons(_vcf_star(tmp_path, samples, rows), max_missing_frac=0.5)
    got = df.set_index("sample")["n_singleton"].to_dict()
    assert got["s1"] == 1 and got["s3"] == 1     # T and G, both private
    assert got["s2"] == 0                        # the deletion
