"""vcf_to_bed: VCF POS is 1-based, BED is not, and that is the entire job.

The test that matters is the round trip -- feed the BED back to `bcftools -R` and it must
select exactly the records it was made from. An off-by-one passes every column-shape check
and fails that one.
"""

import shutil
import subprocess

import pytest

from plasgenomicsutils.lib import vcf_filters as F

pytestmark = pytest.mark.skipif(not shutil.which("bcftools"),
                                reason="bcftools not on PATH")


def _vcf(tmp_path):
    """A SNP, a 4-base deletion, an MNP and another SNP."""
    rows = [(1000, "A", "T"), (2000, "ATTT", "A"), (3000, "AC", "GT"), (4000, "G", "C")]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"]
    for pos, ref, alt in rows:
        hdr.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t222\t.\t.\tGT\t0/1")
    p = tmp_path / "in.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def _rows(path):
    return [ln.split("\t") for ln in open(path).read().splitlines() if ln]


def test_positions_are_zero_based_and_span_the_ref_allele(tmp_path):
    out = str(tmp_path / "o.bed")
    F.vcf_to_bed(_vcf(tmp_path), out)
    r = _rows(out)
    assert [x[1] for x in r] == ["999", "1999", "2999", "3999"]     # POS - 1
    # a SNP is one base; the deletion spans its 4-base REF; the MNP its 2
    assert [int(x[2]) - int(x[1]) for x in r] == [1, 4, 2, 1]


def test_the_bed_selects_exactly_the_records_it_came_from(tmp_path):
    """The round trip. A BED that is one base out still looks perfectly well formed."""
    plain = _vcf(tmp_path)
    vcf = str(tmp_path / "in.vcf.gz")          # -R reads through an index, so make one
    subprocess.run(["bcftools", "view", plain, "-Oz", "-o", vcf], check=True,
                   stderr=subprocess.DEVNULL)
    subprocess.run(["bcftools", "index", "-f", vcf], check=True, stderr=subprocess.DEVNULL)
    bed = str(tmp_path / "o.bed")
    F.vcf_to_bed(vcf, bed, name_column=False)

    def positions(*args):
        out = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\n", *args, vcf],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        return sorted(out.stdout.split("\n"))

    assert positions("-R", bed) == positions()

    # the same check against a deliberately shifted BED, so the assertion above is known to
    # discriminate rather than to pass whatever it is given
    shifted = tmp_path / "shifted.bed"
    shifted.write_text("".join(f"{c}\t{int(s) + 1}\t{int(e) + 1}\n" for c, s, e in _rows(bed)))
    assert positions("-R", str(shifted)) != positions()


def test_snps_only_keeps_substitutions(tmp_path):
    out = str(tmp_path / "o.bed")
    F.vcf_to_bed(_vcf(tmp_path), out, snps_only=True)
    assert [x[1] for x in _rows(out)] == ["999", "3999"]


def test_the_name_column_is_the_canonical_snp_label(tmp_path):
    from plasgenomicsutils.lib.intervals import snp_label

    out = str(tmp_path / "o.bed")
    F.vcf_to_bed(_vcf(tmp_path), out)
    r = _rows(out)
    assert [x[3] for x in r][0] == snp_label("chr1", 999)
    # ...and can be left off for a bare 3-column BED
    bare = str(tmp_path / "b.bed")
    F.vcf_to_bed(_vcf(tmp_path), bare, name_column=False)
    assert all(len(x) == 3 for x in _rows(bare))


def test_no_out_path_writes_to_stdout(tmp_path, capfd):
    F.vcf_to_bed(_vcf(tmp_path))
    out = capfd.readouterr().out
    assert out.startswith("chr1\t999\t1000\t")
    assert len(out.strip().splitlines()) == 4


def test_snp_bed_is_this_restricted_to_snps(tmp_path):
    """The IBD panel writer delegates here, so the two cannot drift apart."""
    a, b = str(tmp_path / "a.bed"), str(tmp_path / "b.bed")
    F.snp_bed(_vcf(tmp_path), a)
    F.vcf_to_bed(_vcf(tmp_path), b, snps_only=True)
    assert open(a).read() == open(b).read()


def _mixed_vcf(tmp_path):
    """Records whose alleles disagree about what type the record is.

    ``bcftools view -v snps`` keeps a record if **any** allele is a substitution, which is
    the trap ``NON_SNP_TYPES`` exists to avoid -- so these are exactly the records that
    tell the two definitions apart.
    """
    rows = [
        (1000, "A", "T"),        # plain SNP                        -> keep
        (2000, "A", "T,ATT"),    # SNP + insertion in one record    -> drop
        (3000, "A", "*,T"),      # SNP under a deletion             -> drop
        (4000, "A", "*"),        # nothing but a spanning deletion  -> drop
        (5000, "AC", "GT"),      # MNP                              -> drop
    ]
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"]
    for pos, ref, alt in rows:
        hdr.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t222\t.\t.\tGT\t0/1")
    p = tmp_path / "mixed.vcf"
    p.write_text("\n".join(hdr) + "\n")
    return str(p)


def test_snps_only_drops_a_record_whose_alleles_are_not_all_substitutions(tmp_path):
    """`A>T,ATT` is not a SNP record, and neither is `A>*,T`; the panel must not say they are.

    ``snp_bed`` runs at the end of every pipeline and produces the IBD SNP panel, so a
    record the callset filter excluded but the panel includes makes the two disagree about
    which positions are SNPs.
    """
    out = str(tmp_path / "o.bed")
    F.vcf_to_bed(_mixed_vcf(tmp_path), out, snps_only=True)
    assert [x[1] for x in _rows(out)] == ["999"]


def test_snps_only_agrees_with_the_callset_filter_on_the_same_records(tmp_path):
    """The panel and the filtered callset must name the same positions."""
    src = _mixed_vcf(tmp_path)
    kept = str(tmp_path / "kept.vcf")
    F.biallelic_snp_filter(src, kept, trim=False, snps_only=True, biallelic=False,
                           mnp_handling="remove")
    expected = subprocess.run(["bcftools", "query", "-f", "%POS0\n", kept],
                              stdout=subprocess.PIPE, text=True,
                              stderr=subprocess.DEVNULL).stdout.split()
    bed = str(tmp_path / "o.bed")
    F.snp_bed(src, bed)
    assert [x[1] for x in _rows(bed)] == expected


def test_a_padded_substitution_is_one_base_at_the_base_that_varies(tmp_path):
    """REF=ATTTA ALT=ATTCA is a SNP -- to the type test, to classify_record, and to the SNP
    filter, which atomises it to T>C at POS+3. The BED used to describe it as a five-base
    interval named for the record's first base, three bases from the one that varies, so the
    panel and the filtered callset disagreed about where the SNP was."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
           "chr1\t1000\t.\tATTTA\tATTCA\t222\t.\t.\tGT\t1/1",
           "chr1\t2000\t.\tG\tC\t222\t.\t.\tGT\t1/1"]
    p = tmp_path / "pad.vcf"
    p.write_text("\n".join(hdr) + "\n")
    bed = tmp_path / "pad.bed"
    F.vcf_to_bed(str(p), str(bed), snps_only=True)
    rows = [ln.split("\t") for ln in bed.read_text().splitlines()]
    assert [(r[0], int(r[1]), int(r[2]), r[3]) for r in rows] == [
        ("chr1", 1002, 1003, "chr1:1002"),     # 1-based 1003: the T>C, one base wide
        ("chr1", 1999, 2000, "chr1:1999"),
    ]
