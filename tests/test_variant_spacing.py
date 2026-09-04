"""variant_spacing: the gaps between consecutive variants, per chromosome.

The two things worth guarding are that a gap never spans a chromosome boundary, and that the
density denominator is the span actually asked about -- including when the requested regions
overlap each other.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib.variant_spacing import (merge_intervals, parse_locus,
                                                   variant_spacing)

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"),
                                reason="bcftools not on PATH")


def _vcf(tmp_path, rows=None, name="in.vcf"):
    """chr1 with variants at 1000/1100/1400, chr2 at 500/2500."""
    rows = rows or [("chr1", 1000), ("chr1", 1100), ("chr1", 1400),
                    ("chr2", 500), ("chr2", 2500)]
    hdr = ["##fileformat=VCFv4.2",
           "##contig=<ID=chr1,length=10000>", "##contig=<ID=chr2,length=20000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"]
    for chrom, pos in rows:
        hdr.append(f"{chrom}\t{pos}\t.\tA\tT\t222\t.\t.\tGT\t0/1")
    p = tmp_path / name
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _by_chrom(rows):
    return {r["chrom"]: r for r in rows}


def test_gaps_never_cross_a_chromosome_boundary(tmp_path):
    """The one mistake this command exists to stop you making by hand."""
    r = _by_chrom(variant_spacing(_vcf(tmp_path), quiet=True))
    assert r["chr1"]["n_variants"] == 3 and r["chr1"]["n_gaps"] == 2   # 100, 300
    assert r["chr2"]["n_variants"] == 2 and r["chr2"]["n_gaps"] == 1   # 2000
    assert r["chr1"]["min"] == 100 and r["chr1"]["max"] == 300
    assert r["chr2"]["min"] == r["chr2"]["max"] == 2000
    # n_gaps totals to variants minus one per chromosome -- never variants minus one
    assert r["all"]["n_gaps"] == r["all"]["n_variants"] - 2


def test_the_all_row_pools_the_gaps_rather_than_averaging_the_rows(tmp_path):
    r = _by_chrom(variant_spacing(_vcf(tmp_path), quiet=True))
    assert r["all"]["min"] == 100          # from chr1
    assert r["all"]["max"] == 2000         # from chr2
    assert r["all"]["p50"] == 300          # median of [100, 300, 2000]


def test_density_uses_the_contig_length_when_no_region_is_given(tmp_path):
    r = _by_chrom(variant_spacing(_vcf(tmp_path), bp_per_cm=1000.0, quiet=True))
    assert r["chr1"]["span_bp"] == 10000
    assert r["chr1"]["variants_per_cm"] == pytest.approx(3 / (10000 / 1000))
    assert r["all"]["span_bp"] == 30000


def test_a_settable_map_rate_moves_only_the_density(tmp_path):
    v = _vcf(tmp_path)
    a = _by_chrom(variant_spacing(v, bp_per_cm=1000.0, quiet=True))["chr1"]
    b = _by_chrom(variant_spacing(v, bp_per_cm=2000.0, quiet=True))["chr1"]
    assert b["variants_per_cm"] == pytest.approx(a["variants_per_cm"] * 2)
    assert a["p50"] == b["p50"]


def _indexed(tmp_path, plain):
    gz = str(tmp_path / "in.vcf.gz")
    subprocess.run(["bcftools", "view", plain, "-Oz", "-o", gz], check=True,
                   stderr=subprocess.DEVNULL)
    subprocess.run(["bcftools", "index", "-f", gz], check=True, stderr=subprocess.DEVNULL)
    return gz


def test_overlapping_regions_are_not_counted_twice_in_the_denominator(tmp_path):
    """bcftools reports each variant once however the regions overlap, so the span has to
    agree with it -- otherwise a nested region silently deflates the density."""
    v = _indexed(tmp_path, _vcf(tmp_path))
    one = _by_chrom(variant_spacing(v, locus="chr1:1000-1400", quiet=True))["chr1"]
    two = _by_chrom(variant_spacing(v, locus="chr1:1000-1400,chr1:1050-1200",
                                    quiet=True))["chr1"]
    assert one["span_bp"] == 401                      # 1-based inclusive, both ends counted
    assert two["span_bp"] == one["span_bp"]
    assert two["n_variants"] == one["n_variants"] == 3


def test_a_locus_is_one_based_and_a_bed_is_zero_based(tmp_path):
    """The two conventions bcftools itself uses, kept rather than unified -- so each reads
    the way its own format does. The same interval written both ways must agree."""
    v = _indexed(tmp_path, _vcf(tmp_path))
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t999\t1400\n")               # 0-based half-open == chr1:1000-1400
    a = _by_chrom(variant_spacing(v, locus="chr1:1000-1400", quiet=True))["chr1"]
    b = _by_chrom(variant_spacing(v, bed=str(bed), quiet=True))["chr1"]
    assert a["span_bp"] == b["span_bp"] and a["n_variants"] == b["n_variants"]


def test_a_region_with_no_variants_is_an_answer_not_an_error(tmp_path):
    v = _indexed(tmp_path, _vcf(tmp_path))
    r = _by_chrom(variant_spacing(v, locus="chr1:5000-6000", quiet=True))["chr1"]
    assert r["n_variants"] == 0 and r["n_gaps"] == 0
    assert r["span_bp"] == 1001 and r["variants_per_cm"] == 0.0
    assert r["p50"] is None            # undefined, not zero


def test_one_variant_on_a_chromosome_has_no_gap_statistics(tmp_path):
    v = _vcf(tmp_path, rows=[("chr1", 1000), ("chr2", 500), ("chr2", 900)])
    r = _by_chrom(variant_spacing(v, quiet=True))
    assert r["chr1"]["n_variants"] == 1 and r["chr1"]["n_gaps"] == 0
    assert r["chr1"]["p50"] is None and r["chr1"]["min"] is None


def test_an_unknown_contig_is_rejected_by_name(tmp_path):
    v = _indexed(tmp_path, _vcf(tmp_path))
    with pytest.raises(SystemExit, match="chr9 is not a contig"):
        variant_spacing(v, locus="chr9:1-100", quiet=True)


def test_locus_parsing(tmp_path):
    assert parse_locus("c:100-200") == [("c", 99, 200)]
    assert parse_locus("c") == [("c", 0, -1)]                   # bare contig
    # the comma separates regions, so it cannot group digits too -- and the error says so
    with pytest.raises(SystemExit, match="cannot contain commas"):
        parse_locus("c:1,000-2,000")
    assert merge_intervals(parse_locus("c:100-200,c:150-160")) == {"c": [(99, 200)]}
    with pytest.raises(SystemExit, match="ends before it starts"):
        parse_locus("c:200-100")


def test_locus_and_bed_together_is_refused(tmp_path):
    with pytest.raises(SystemExit, match="not both"):
        variant_spacing(_vcf(tmp_path), locus="chr1", bed="x.bed", quiet=True)
