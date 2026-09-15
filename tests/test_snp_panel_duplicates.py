"""A SNP panel names positions, so one position must appear once.

`bcftools norm -m-` splits a multiallelic record into one record per ALT **at the same
position**. Fed to a panel loader that appends one row per line, that yields several panel
entries for one SNP, and every count built on the panel -- the panel length, the IBD matrix
columns a block covers, the variant density per cM -- is inflated by however many alternates
the callset happened to carry.

The merged record is this package's interchange form (MULTIALLELIC_PLAN.md, decision 1), so
split input is a mistake rather than a supported shape. It should say so: silently
de-duplicating would let the rest of the run proceed on a panel that does not match the
callset it came from.
"""

import pytest

from plasgenomicsutils.lib.vcf_io import SnpPanel, positions_frame


def _split_vcf(tmp_path):
    """What `bcftools norm -m-` makes of one triallelic record."""
    p = tmp_path / "split.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1000\t.\tA\tC\t.\t.\t.\n"
        "chr1\t1000\t.\tA\tG\t.\t.\t.\n"
        "chr1\t2000\t.\tA\tT\t.\t.\t.\n")
    return str(p)


def _merged_vcf(tmp_path):
    p = tmp_path / "merged.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1000\t.\tA\tC,G\t.\t.\t.\n"
        "chr1\t2000\t.\tA\tT\t.\t.\t.\n")
    return str(p)


def test_a_split_callset_is_refused_by_the_panel_loader(tmp_path):
    with pytest.raises(SystemExit, match="more than once"):
        SnpPanel.from_vcf(_split_vcf(tmp_path))


def test_the_message_names_the_position_and_the_remedy(tmp_path):
    with pytest.raises(SystemExit) as e:
        SnpPanel.from_vcf(_split_vcf(tmp_path))
    msg = str(e.value)
    assert "chr1:999" in msg
    assert "norm -m" in msg, "it should say how to put the callset back together"


def test_the_merged_form_loads_and_counts_one_snp_per_position(tmp_path):
    panel = SnpPanel.from_vcf(_merged_vcf(tmp_path))
    assert len(panel) == 2
    assert panel.labels == ["chr1:999", "chr1:1999"]
    assert list(panel.snps_in_block("chr1", 0, 5000)) == [0, 1]


def test_a_bed_panel_is_checked_the_same_way(tmp_path):
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t999\t1000\nchr1\t999\t1000\nchr1\t1999\t2000\n")
    with pytest.raises(SystemExit, match="more than once"):
        SnpPanel.from_bed(str(bed))


def test_positions_frame_refuses_split_input_too(tmp_path):
    """It feeds variant density, where a duplicate position reads as a zero-length gap."""
    with pytest.raises(SystemExit, match="more than once"):
        positions_frame(_split_vcf(tmp_path), "vcf")
    df = positions_frame(_merged_vcf(tmp_path), "vcf")
    assert len(df) == 2
