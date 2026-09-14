"""ld_recombination's own logic -- the haploid-panel prep and the pyrho input encoding.

pyrho itself is not exercised here (it lives in a separate env); these pin the parts we
own: rescaling read_dosages' 0/2/-1 to a complete 0/1 panel, and the pseudo-diploid VCF
encoding that pyrho's --ploidy 1 flattens back into independent haplotypes."""

from __future__ import annotations

import numpy as np
import pytest

from plasgenomicsutils.lib.fws import monoclonal_samples
from plasgenomicsutils.scripts.ld.recombination import (
    build_haplotype_panel,
    parse_ldhat_res,
    read_pyrho_map,
    resolve_ldhat,
    resolve_pyrho_cmd,
    write_pyrho_vcf,
)


def test_panel_rescales_dosage_and_splits_by_chromosome():
    # two chromosomes; dosage 0/2 -> 0/1, all sites polymorphic and complete
    gn = np.array([[0, 2, 2, 0], [2, 2, 0, 0], [0, 0, 2, 2]], dtype=np.int8)
    chrom = np.array(["c1", "c1", "c2"])
    pos = np.array([10, 20, 30])
    out = dict((c, (h, p)) for c, h, p in build_haplotype_panel(
        gn, chrom, pos, max_missing=0.0, maf=0.0))
    assert set(out) == {"c1", "c2"}
    assert set(np.unique(out["c1"][0])) <= {0, 1}          # rescaled to 0/1
    assert out["c1"][0].shape == (2, 4) and out["c2"][0].shape == (1, 4)


def test_monomorphic_and_over_missing_sites_are_dropped():
    # row0 monomorphic (all alt), row1 half missing, row2 good
    gn = np.array([[2, 2, 2, 2], [0, 2, -1, -1], [0, 0, 2, 2]], dtype=np.int8)
    chrom = np.array(["c1", "c1", "c1"])
    pos = np.array([10, 20, 30])
    (c, hap, p), = list(build_haplotype_panel(gn, chrom, pos, max_missing=0.10, maf=0.0))
    assert list(p) == [30]                                  # only the clean polymorphic site


def test_missing_calls_are_imputed_to_the_major_allele():
    # one site, 4 samples: three ref (0) one missing -> imputes to 0, then monomorphic->dropped
    gn = np.array([[0, 0, 0, -1]], dtype=np.int8)
    assert list(build_haplotype_panel(gn, np.array(["c"]), np.array([1]),
                                      max_missing=0.5, maf=0.0)) == []
    # now with a real minor allele: missing imputes to major (0), site stays polymorphic
    gn = np.array([[0, 0, 2, -1]], dtype=np.int8)
    (_, hap, _), = list(build_haplotype_panel(gn, np.array(["c"]), np.array([1]),
                                              max_missing=0.5, maf=0.0))
    assert list(hap[0]) == [0, 0, 1, 0]


def test_pyrho_vcf_pairs_haplotypes_and_flatten_is_lossless(tmp_path):
    # 4 haplotypes over 2 SNPs; pairing then flattening a|b must reproduce the columns
    hap = np.array([[0, 1, 1, 0], [1, 1, 0, 0]], dtype=np.int8)
    pos = np.array([100, 200])
    out = tmp_path / "in.vcf"
    n_hap, odd = write_pyrho_vcf("chrX", hap, pos, str(out))
    assert (n_hap, odd) == (4, 0)
    lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
    # reconstruct haplotypes by flattening the a|b genotypes, per SNP row
    for j, line in enumerate(lines):
        gts = line.split("\t")[9:]
        flat = [int(x) for g in gts for x in g.split("|")]
        assert flat == list(hap[j])
    assert "##contig=<ID=chrX>" in out.read_text()          # silences pyrho's warning
    assert lines[0].split("\t")[1] == "101"                 # 0-based pos -> 1-based VCF


def test_odd_panel_drops_one_haplotype_for_pairing(tmp_path):
    hap = np.array([[0, 1, 1]], dtype=np.int8)              # 3 haplotypes
    out = tmp_path / "odd.vcf"
    n_hap, odd = write_pyrho_vcf("c", hap, np.array([5]), str(out))
    assert (n_hap, odd) == (2, 1)
    gts = [l for l in out.read_text().splitlines() if not l.startswith("#")][0].split("\t")[9:]
    assert len(gts) == 1 and gts[0] == "0|1"                # one pair, last hap dropped


def test_monoclonal_samples_thresholds_on_fws(tmp_path):
    t = tmp_path / "fws.tsv"
    t.write_text("sample\tfws\tn_sites\tmonoclonal\n"
                 "s1\t0.99\t100\tTrue\n"
                 "s2\t0.80\t100\tFalse\n"
                 "s3\t0.96\t100\tTrue\n")
    assert monoclonal_samples(str(t), fws_min=0.95) == ["s1", "s3"]
    assert monoclonal_samples(str(t), fws_min=0.70) == ["s1", "s2", "s3"]


def test_monoclonal_samples_falls_back_to_the_boolean_column(tmp_path):
    t = tmp_path / "mono.tsv"                                # no fws column
    t.write_text("sample\tmonoclonal\ns1\tTrue\ns2\tfalse\ns3\t1\n")
    assert monoclonal_samples(str(t)) == ["s1", "s3"]


def test_read_pyrho_map_parses_headerless_rows(tmp_path):
    t = tmp_path / "out.tsv"
    t.write_text("100\t200\t0.0011\n200\t300\t0.0022\n")
    df = read_pyrho_map(str(t), "c7")
    assert list(df.columns) == ["chrom", "start", "end", "rho_per_bp"]
    assert list(df["chrom"]) == ["c7", "c7"] and df["rho_per_bp"].iloc[1] == pytest.approx(0.0022)


def test_resolve_pyrho_cmd_explicit_and_auto():
    assert resolve_pyrho_cmd("mamba run -n foo pyrho") == ["mamba", "run", "-n", "foo", "pyrho"]
    auto = resolve_pyrho_cmd("auto")                        # PATH pyrho or a manager fallback
    assert auto[-1] == "pyrho"


def test_resolve_ldhat_finds_binaries_and_handles_arch_prefix(tmp_path):
    for name in ("interval", "lkgen", "stat"):
        (tmp_path / name).write_text("")                    # stub binaries
    prefix, bins = resolve_ldhat(str(tmp_path), arch_prefix="none")
    assert prefix == [] and bins["interval"].endswith("interval")
    prefix, _ = resolve_ldhat(str(tmp_path), arch_prefix="arch -x86_64")
    assert prefix == ["arch", "-x86_64"]


def test_resolve_ldhat_errors_without_dir_or_binaries(tmp_path):
    with pytest.raises(SystemExit):
        resolve_ldhat(None)
    with pytest.raises(SystemExit):                         # dir exists but no binaries
        resolve_ldhat(str(tmp_path))


def test_parse_ldhat_res_drops_summary_row(tmp_path):
    # LDhat stat res.txt: header, then a locus -1 total row, then per-locus rate rows
    t = tmp_path / "res.txt"
    t.write_text("Loci\tMean_rho\tMedian\tL95\tU95\n"
                 "-1.000\t40098.0\t40000\t38000\t42000\n"
                 "77.130\t55.3\t56.0\t52\t60\n"
                 "80.500\t12.4\t13.0\t10\t15\n")
    kb, rho = parse_ldhat_res(str(t))
    assert list(kb) == [77.130, 80.500]                     # the -1 summary row is gone
    assert list(rho) == [55.3, 12.4]
