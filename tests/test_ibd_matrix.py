"""Known-answer + invariant tests for the IBD matrix layer.

Self-contained: builds tiny inputs in a tmp dir. No external data. When the
public test dataset lands, add golden-file fixtures alongside these.
"""

import numpy as np
import pandas as pd

from plasgenomicsutils.lib import ibd_matrix
from plasgenomicsutils.lib.vcf_io import SnpPanel


def _write(path, text):
    path.write_text(text)
    return str(path)


def test_load_matrix_accepts_prefix_or_npz(tmp_path):
    from scipy.sparse import csr_matrix
    prefix = str(tmp_path / "ibd_matrix")
    ibd_matrix.save_matrix(csr_matrix(np.array([[1, 0], [1, 1]])),
                           ["p1__p2", "p3__p4"], ["c:1", "c:2"], prefix)
    a = ibd_matrix.load_matrix(prefix)              # the prefix, as documented
    b = ibd_matrix.load_matrix(prefix + ".npz")     # the full .npz path
    assert (a[0] != b[0]).nnz == 0 and a[1] == b[1] == ["p1__p2", "p3__p4"]


def test_build_matrix_known_answer(tmp_path):
    # 4 SNPs on chr1 at 1-based POS 100,200,300,400 -> pos0 99,199,299,399
    vcf = _write(tmp_path / "s.vcf",
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tT\t.\tPASS\t.\n"
        "chr1\t200\t.\tA\tT\t.\tPASS\t.\n"
        "chr1\t300\t.\tA\tT\t.\tPASS\t.\n"
        "chr1\t400\t.\tA\tT\t.\tPASS\t.\n")
    # one IBD block for pair (A,B) covering pos0 [150,350] -> SNPs at 199,299 (cols 1,2)
    blocks = _write(tmp_path / "b.tsv",
        "sample1\tsample2\tchr\tstart\tend\tdifferent\tNsnp\n"
        "A\tB\tchr1\t150\t350\t0\t2\n"
        "A\tB\tchr1\t150\t250\t0\t1\n"      # overlapping second block -> still binary
        "B\tC\tchr1\t0\t1000\t1\t9\n")       # different==1 -> excluded, but pair exists

    panel = SnpPanel.load(vcf, "vcf")
    assert panel.labels == ["chr1:99", "chr1:199", "chr1:299", "chr1:399"]

    blocks_df = ibd_matrix.read_blocks(blocks)
    pair_to_row, pair_labels = ibd_matrix.build_pair_index(blocks_df)
    # pairs come from ALL blocks (incl. different==1): A__B and B__C
    assert pair_labels == ["A__B", "B__C"]

    mat = ibd_matrix.build_matrix(blocks_df, panel, pair_to_row).toarray()
    # A__B row: cols 1,2 set (from the different==0 blocks), binary despite overlap
    assert mat[pair_to_row[("A", "B")]].tolist() == [0, 1, 1, 0]
    # B__C only had a different==1 block -> excluded -> all zeros, but row exists
    assert mat[pair_to_row[("B", "C")]].tolist() == [0, 0, 0, 0]
    assert set(np.unique(mat)).issubset({0, 1})  # strictly binary


def test_snps_in_block_binary_search(tmp_path):
    bed = _write(tmp_path / "s.bed",
        "chr1\t10\t11\nchr1\t20\t21\nchr1\t30\t31\nchr2\t10\t11\n")
    panel = SnpPanel.load(bed, "bed")
    # half-open [20, 30) on chr1 -> only the SNP at 20; 30 is the excluded end
    hits = sorted(panel.snps_in_block("chr1", 20, 30).tolist())
    assert sorted(panel.df.loc[hits, "pos0"].tolist()) == [20]
    # extending the end past 30 takes it in
    hits = sorted(panel.snps_in_block("chr1", 20, 31).tolist())
    assert sorted(panel.df.loc[hits, "pos0"].tolist()) == [20, 30]
    # chromosome not present -> empty
    assert panel.snps_in_block("chrX", 0, 100).size == 0


def test_pair_summary_fraction():
    from scipy.sparse import csr_matrix
    mat = csr_matrix(np.array([[1, 0, 1, 0], [0, 0, 0, 0]], dtype=np.uint8))
    summ = ibd_matrix.pair_summary(mat, ["A__B", "A__C"])
    assert summ["n_ibd_snps"].tolist() == [2, 0]
    assert summ["frac_ibd"].tolist() == [0.5, 0.0]
