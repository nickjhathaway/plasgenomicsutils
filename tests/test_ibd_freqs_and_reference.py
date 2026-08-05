"""Tests for single-pass allele frequencies and the reference registry."""

import numpy as np

from plasgenomicsutils.lib.ibd_freqs import compute_allele_freqs
from plasgenomicsutils.lib import reference as R


def test_reference_registry():
    ref = R.get_reference("pf3d7")
    assert ref.species == "Plasmodium falciparum"
    assert len(ref.core_chrom_lengths_bp) == 14
    assert ref.bp_per_cm == 15000.0
    assert "pf3d7" in R.available_references()


def test_normalise_chr():
    assert R.normalise_chr("Pf3D7_07_v3") == "7"
    assert R.normalise_chr("chr07") == "7"
    assert R.normalise_chr("14") == "14"
    assert R.normalise_chr(3) == "3"


def test_unknown_reference_errors():
    import pytest
    with pytest.raises(SystemExit):
        R.get_reference("does_not_exist")


def test_single_pass_allele_freqs(tmp_path):
    # pysam is required for this test; skip cleanly if unavailable
    import pytest
    pytest.importorskip("pysam")

    vcf = tmp_path / "s.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        # S1,S2 region A ; S3,S4 region B
        "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\t0/0\t1/1\t1/1\n"   # global 6/8, A 2/4, B 4/4
        "chr1\t200\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t./.\n")  # global 1/6, A 0/4, B 1/2

    s2g = {"S1": "A", "S2": "A", "S3": "B", "S4": "B"}
    gdf, rdf = compute_allele_freqs(str(vcf), sample_to_group=s2g, zero_based=False)

    g = dict(zip(gdf["snp_id"], gdf["af"]))
    assert np.isclose(g["chr1:100"], 6 / 8)
    assert np.isclose(g["chr1:200"], 1 / 6)

    rA = rdf[(rdf["group"] == "A")].set_index("snp_id")["af"].to_dict()
    rB = rdf[(rdf["group"] == "B")].set_index("snp_id")["af"].to_dict()
    assert np.isclose(rA["chr1:100"], 2 / 4)
    assert np.isclose(rB["chr1:100"], 4 / 4)
    assert np.isclose(rA["chr1:200"], 0 / 4)
    assert np.isclose(rB["chr1:200"], 1 / 2)  # S4 is ./. so only S3 counts (1 alt / 2)

    # group table is group-major, sorted groups
    assert list(rdf["group"].unique()) == ["A", "B"]
