"""Metadata tables are read with the delimiter auto-detected, so a .tsv `--meta`
(the usual case) works without a separator flag — regression for analyze_ibd_matrix /
ibd_selection_statistic reading a tab file as comma-CSV."""

import pandas as pd

from plasgenomicsutils.utils.small_utils import Utils

TAB = "sample\tsite\tregion\tyear\nIM102\tAduku\tAduku\t2023\nIM113\tMetu\tNorthwest\t2023\n"
CSV = "sample,region\nIM102,Aduku\nIM113,Northwest\n"


def test_read_table_autodetects_tab(tmp_path):
    p = tmp_path / "meta.tsv"
    p.write_text(TAB)
    m = Utils.read_table(str(p))
    assert list(m.columns) == ["sample", "site", "region", "year"]
    # the operation that used to crash (KeyError: 'sample') now works
    assert m.set_index("sample")["region"].to_dict() == {"IM102": "Aduku", "IM113": "Northwest"}


def test_read_table_autodetects_comma(tmp_path):
    p = tmp_path / "meta.csv"
    p.write_text(CSV)
    m = Utils.read_table(str(p))
    assert list(m.columns) == ["sample", "region"]


def test_read_table_explicit_sep_override(tmp_path):
    p = tmp_path / "meta.tsv"
    p.write_text(TAB)
    m = Utils.read_table(str(p), sep="\t")
    assert "sample" in m.columns
