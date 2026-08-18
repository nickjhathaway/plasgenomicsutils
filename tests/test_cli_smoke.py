"""End-to-end runs of the commands the rest of the suite exercises only at library level.

Each drives the CLI the way a user does -- argv in, files out -- so a leaf module that
wires an argument to the wrong parameter fails here even when the library function it
calls is well covered.
"""

import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

needs_bcftools = pytest.mark.skipif(shutil.which("bcftools") is None,
                                    reason="bcftools not on PATH")


def _read_maybe_gz(path):
    """The commands compress their tabular output whatever the name says."""
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw).decode()
    return raw.decode()


def run(*argv):
    p = subprocess.run([sys.executable, "-m", "plasgenomicsutils.cli", *map(str, argv)],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert p.returncode == 0, p.stdout
    return p.stdout


def _vcf(path, n_samples=6, n_sites=40, with_ad=True):
    import random
    random.seed(3)
    samples = [f"s{i}" for i in range(1, n_samples + 1)]
    fmt = "GT:AD" if with_ad else "GT"
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=Pf3D7_01_v3,length=640851>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">']
    if with_ad:
        hdr.append('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">')
        hdr.append('##FORMAT=<ID=ADS,Number=1,Type=Integer,Description="ADS">')
        fmt = "GT:AD:ADS"
    hdr.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
    rows = []
    for i in range(n_sites):
        cells = []
        for _ in samples:
            alt = random.random() < 0.35
            d = random.randint(20, 60)
            a = d if alt else 0
            g = "1/1" if alt else "0/0"
            cells.append(f"{g}:{d - a},{a}:{d}" if with_ad else g)
        rows.append(f"Pf3D7_01_v3\t{1000 * (i + 1)}\t.\tA\tT\t100\t.\t.\t{fmt}\t"
                    + "\t".join(cells))
    Path(path).write_text("\n".join(hdr + rows) + "\n")
    return str(path)


@needs_bcftools
def test_calculate_fws_writes_one_row_per_sample(tmp_path):
    vcf = _vcf(tmp_path / "in.vcf")
    out = tmp_path / "fws.tsv"
    run("calculate_fws", "--input-vcf", vcf, "--out", out, "--monoclonal-threshold", 0.95)

    rows = [l.split("\t") for l in out.read_text().strip().splitlines()]
    assert rows[0][:2] == ["sample", "fws"]
    assert len(rows) - 1 == 6
    fws = [float(r[1]) for r in rows[1:]]
    assert all(0 <= f <= 1 for f in fws)          # Fws is a proportion
    # the monoclonal flag has to agree with the threshold it was given
    mono = {r[0]: r[-1] for r in rows[1:]}
    for r in rows[1:]:
        assert (float(r[1]) >= 0.95) == (mono[r[0]].lower() in ("true", "1"))


@needs_bcftools
def test_locus_missingness_filter_drops_thinly_covered_loci(tmp_path):
    vcf = _vcf(tmp_path / "in.vcf")
    keep = tmp_path / "keep.bcf"
    drop = tmp_path / "drop.bcf"
    run("locus_missingness_filter", "--input", vcf, "--output", keep, "--ads-min", 1)
    run("locus_missingness_filter", "--input", vcf, "--output", drop, "--ads-min", 999)

    def n(p):
        return len(subprocess.run(["bcftools", "view", "-H", str(p)],
                                  stdout=subprocess.PIPE, text=True).stdout.splitlines())
    assert n(keep) == 40          # every locus is deeply covered in the fixture
    assert n(drop) == 0           # nothing clears an impossible ADS floor


def test_ibd_fraction_and_snp_density_writes_both_tables(tmp_path):
    blocks = tmp_path / "blocks.txt"
    blocks.write_text(
        "sample1\tsample2\tchr\tstart\tend\tdifferent\tNsnp\n"
        "s1\ts2\t1\t100000\t200000\t0\t50\n"
        "s1\ts3\t1\t300000\t340000\t0\t30\n")
    snps = tmp_path / "snps.bed"
    snps.write_text("".join(f"Pf3D7_01_v3\t{p}\t{p + 1}\n" for p in range(100000, 400000, 1000)))
    out = tmp_path / "frac"
    run("ibd_fraction_and_snp_density", "--blocks", blocks, "--snps", snps,
        "--snp-format", "bed", "--output", out)

    pair = tmp_path / "frac.pair_ibd_fraction.tsv.gz"
    assert pair.exists()
    with gzip.open(pair, "rt") as fh:
        head = fh.readline().split("\t")
    for col in ("sample1", "sample2", "pair", "ibd_fraction_accessible"):
        assert col in [h.strip() for h in head], f"{col} missing from {head}"


def test_coverage_dropout_regions_finds_the_window_empty_in_everyone(tmp_path):
    win = tmp_path / "windows.tsv"
    lines = ["sample\tchrom\tstart\tend\tmean_depth"]
    for s in range(1, 6):
        for i, start in enumerate(range(0, 400000, 100000)):
            # the third window is empty in every sample; the rest are well covered
            depth = 0.5 if i == 2 else 40
            lines.append(f"s{s}\tPf3D7_01_v3\t{start}\t{start + 100000}\t{depth}")
    win.write_text("\n".join(lines) + "\n")

    out = tmp_path / "dropout.tsv"
    bed = tmp_path / "dropout.bed"
    run("coverage_dropout_regions", "--windows", win, "--output", out,
        "--bed-output", bed, "--min-depth", 5, "--min-frac-samples", 0.9)

    text = _read_maybe_gz(out)
    body = [l for l in text.strip().splitlines()[1:] if l]
    assert len(body) == 1, text
    fields = body[0].split("\t")
    assert "200000" in fields and "300000" in fields
    assert _read_maybe_gz(bed).strip().split("\t")[:3] == ["Pf3D7_01_v3", "200000", "300000"]


def test_strand_bias_scan_flags_a_one_sided_alt(tmp_path):
    """A real variant has alt reads on both strands; an SSE artifact has them on one."""
    hdr = ["##fileformat=VCFv4.2", "##contig=<ID=Pf3D7_01_v3,length=640851>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           '##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="fwd">',
           '##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="rev">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"]
    rows = [
        # balanced: 20 alt reads each way in both samples
        "Pf3D7_01_v3\t1000\t.\tA\tT\t100\t.\t.\tGT:ADF:ADR\t0/1:20,20:20,20\t0/1:20,20:20,20",
        # one-sided: every alt read on the forward strand
        "Pf3D7_01_v3\t2000\t.\tA\tT\t100\t.\t.\tGT:ADF:ADR\t0/1:20,40:20,0\t0/1:20,40:20,0",
    ]
    vcf = tmp_path / "sb.vcf"
    vcf.write_text("\n".join(hdr + rows) + "\n")

    tsv, bed = tmp_path / "sb.tsv", tmp_path / "sb.bed"
    run("strand_bias_scan", "--input-vcf", vcf, "--out-tsv", tsv, "--out-bed", bed,
        "--min-drop-samples", 2)

    lines = _read_maybe_gz(tsv).strip().splitlines()
    assert len(lines) - 1 == 4                      # 2 sites x 2 samples
    # only the one-sided position reaches the blacklist
    bed_txt = _read_maybe_gz(bed).strip()
    assert "1999" in bed_txt or "2000" in bed_txt
    assert bed_txt.count("\n") == 0                 # exactly one region


def _tiny_bam(path, contig="Pf3D7_01_v3", length=10000, n_reads=200, read_len=100):
    """A BAM of uniformly tiled 100 bp reads, so depth is predictable."""
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": contig, "LN": length}],
              "RG": [{"ID": "rg1", "SM": "s1"}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(n_reads):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * read_len
            a.flag = 0
            a.reference_id = 0
            a.reference_start = (i * (length // n_reads)) % (length - read_len)
            a.mapping_quality = 60
            a.cigarstring = f"{read_len}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            a.set_tag("RG", "rg1")
            out.write(a)
    pysam.sort("-o", str(path) + ".sorted", str(path))
    Path(str(path) + ".sorted").replace(path)
    pysam.index(str(path))
    return str(path)


def test_coverage_depth_stats_writes_per_chromosome_and_genome_rows(tmp_path):
    bam = _tiny_bam(tmp_path / "s1.bam")
    out = tmp_path / "cov.tsv"
    win = tmp_path / "cov_windows.tsv"
    run("coverage_depth_stats", "--bam", bam, "--thresholds", 1, "--window", 1000,
        "--engine", "pysam", "--output", out, "--windows-output", win)

    lines = _read_maybe_gz(out).strip().splitlines()
    head = lines[0].split("\t")
    for col in ("sample", "chrom", "mean", "engine"):
        assert col in head, head
    rows = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    assert any(r["chrom"] == "genome" for r in rows), "no genome-wide row"
    assert {r["engine"] for r in rows} == {"pysam"}          # the engine is recorded
    assert all(float(r["mean"]) > 0 for r in rows)

    wlines = _read_maybe_gz(win).strip().splitlines()
    assert len(wlines) > 1
    for col in ("start", "end", "mean_depth"):
        assert col in wlines[0].split("\t")


def test_strand_read_check_tabulates_the_reads_over_one_site(tmp_path):
    """The per-read table behind a strand-bias call: one row per read at the position."""
    bam = _tiny_bam(tmp_path / "s1.bam")
    prefix = tmp_path / "site"
    run("strand_read_check", "--bam", bam, "--pos", "Pf3D7_01_v3:1050",
        "--ref-base", "A", "--out", prefix)

    written = sorted(p.name for p in tmp_path.glob("site*"))
    assert written, f"nothing written; tmp holds {[p.name for p in tmp_path.iterdir()]}"
    table = next(p for p in tmp_path.glob("site*") if p.suffix in (".tsv", ".txt", ".gz"))
    head = _read_maybe_gz(table).splitlines()[0].split("\t")
    # the columns the artifact call is read off: which strand, which base, where in the read
    for col in ("strand", "base"):
        assert any(col in h for h in head), head
