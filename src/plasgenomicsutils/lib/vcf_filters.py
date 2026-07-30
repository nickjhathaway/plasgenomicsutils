"""Parameterized VCF/BCF filtering steps, each backed by bcftools/bedtools.

Every step takes an input and output path plus tunable thresholds, so the same
building blocks can be run individually or chained by ``filter_pipeline``.
"""

from __future__ import annotations

import subprocess

from .bcftools import out_flag, q, require, sh


def hard_qc_filter(inp: str, out: str, *, qd: float = 20, mq: float = 55,
                   sor: float = 3, mqranksum: float = -5.0,
                   readposranksum: float = -5.0, fs: float | None = None) -> None:
    """GATK-style hard filter on INFO metrics; keep records still flagged PASS."""
    parts = [
        f"QD < {qd}",
        f"MQ < {mq}",
        f"SOR > {sor}",
        f'(MQRankSum!="." && MQRankSum < {mqranksum})',
        f'(ReadPosRankSum!="." && ReadPosRankSum < {readposranksum})',
    ]
    if fs is not None:
        parts.append(f"FS > {fs}")
    expr = " || ".join(parts)
    fmt = out_flag(out)
    sh(f"bcftools filter -m + -s FAIL -e {q(expr)} {q(inp)} -Ou "
       f"| bcftools view -f PASS -O{fmt} -o {q(out)}", tools=("bcftools",))


def singleton_add_ads(inp: str, out: str, *, min_samples: int = 1) -> None:
    """Drop variants seen as ALT in <= min_samples samples; add FORMAT/ADS.

    ADS (summed allelic depth, the reads actually used to genotype) is a more
    honest per-site depth than DP. It is the sum over *all* AD entries via
    ``smpl_sum``, so it stays correct on multiallelic sites — the older
    ``AD[*:0]+AD[*:1]`` form counted only ref + first ALT and silently
    undercounted anything with a second ALT allele.
    """
    fmt = out_flag(out)
    # Single '&' is required inside COUNT(): it combines the conditions per
    # sample, so COUNT() tallies samples meeting both. '&&' would collapse them
    # across the whole record and match every site.
    keep = f'COUNT(GT!="RR" & GT!="mis") > {min_samples}'
    ads = "FORMAT/ADS=int(smpl_sum(FORMAT/AD))"
    sh(f"bcftools view -i {q(keep)} {q(inp)} -Ou "
       f"| bcftools +fill-tags -O{fmt} -o {q(out)} -- -t {q(ads)}", tools=("bcftools",))


def biallelic_snp_filter(inp: str, out: str, *, trim: bool = True) -> None:
    """Keep biallelic SNPs, optionally trimming unused ALT alleles first.

    When ``trim`` is set (default), ALT alleles absent from every genotype are
    dropped before the biallelic test (``bcftools view --trim-alt-alleles``). Run
    this AFTER re-genotyping (``filter_ad_regenotype``): a site that looked
    multiallelic only because one sample carried a low-level artifact allele —
    now re-genotyped away — collapses to a genuine biallelic SNP and is kept,
    instead of being discarded by a naive ``-m2 -M2``. Sites left with no ALT
    after trimming are dropped by ``-v snps``.
    """
    fmt = out_flag(out)
    trim_view = "bcftools view --trim-alt-alleles {} -Ou | ".format(q(inp)) if trim else ""
    src = "-" if trim else q(inp)
    sh(f"{trim_view}bcftools view -m2 -M2 -v snps {src} -O{fmt} -o {q(out)}", tools=("bcftools",))


def region_filter(inp: str, out: str, *, bed: str, exclude: bool) -> None:
    """Keep (``exclude=False``) or drop (``exclude=True``) variants overlapping ``bed``.

    Backs the core-genome (keep), tandem-repeat and paralog (drop) masks. bedtools
    emits VCF text which is piped back through bcftools so any output format works.
    """
    fmt = out_flag(out)
    v = "-v " if exclude else ""
    sh(f"bcftools view {q(inp)} "
       f"| bedtools intersect {v}-header -a stdin -b {q(bed)} "
       f"| bcftools view -O{fmt} -o {q(out)}", tools=("bcftools", "bedtools"))


def tandem_repeat_mask(inp: str, out: str, *, bed: str) -> None:
    """Remove variants overlapping a tandem-repeat BED (artifact-prone regions)."""
    region_filter(inp, out, bed=bed, exclude=True)


def paralog_mask(inp: str, out: str, *, bed: str) -> None:
    """Remove variants overlapping paralogous / multigene-family genes (mismapping-prone)."""
    region_filter(inp, out, bed=bed, exclude=True)


def core_region_filter(inp: str, out: str, *, bed: str) -> None:
    """Keep only variants inside the core-genome BED (drop subtelomeric/hypervariable)."""
    region_filter(inp, out, bed=bed, exclude=False)


def sample_coverage_filter(inp: str, out: str, *, ads_min: int = 10,
                           frac_min: float = 0.85,
                           dropped_samples_path: str | None = None) -> list[str]:
    """Drop samples covered (ADS >= ads_min) at < frac_min of loci; refresh AC/AN/AF.

    Returns the list of dropped sample names.
    """
    require("bcftools")
    total = _count(inp)
    counts = _per_sample_covered(inp, ads_min)
    all_samples = _samples(inp)
    dropped = sorted(s for s in all_samples if (counts.get(s, 0) / total if total else 0) < frac_min)

    if dropped_samples_path:
        with open(dropped_samples_path, "w") as fh:
            fh.write("\n".join(dropped) + ("\n" if dropped else ""))

    fmt = out_flag(out)
    # The singleton re-filter always runs: re-genotyping upstream can leave sites
    # supported by a single sample even when no samples are dropped here.
    keep = 'COUNT(GT!="RR" & GT!="mis") > 1'
    drop_arg = (dropped_samples_path or _write_tmp_list(dropped)) if dropped else None
    view_in = (f"bcftools view -S ^{q(drop_arg)} {q(inp)} -Ou"
               if drop_arg else f"bcftools view {q(inp)} -Ou")
    sh(f"{view_in} "
       f"| bcftools view -i {q(keep)} -Ou "
       f"| bcftools +fill-tags -O{fmt} -o {q(out)} -- -t AC,AN,AF", tools=("bcftools",))
    return dropped


def locus_missingness_filter(inp: str, out: str, *, f_missing_max: float = 0.05,
                             ads_min: int = 10, sample_frac_min: float = 0.95) -> None:
    """Keep loci with < f_missing_max missing AND >= sample_frac_min at ADS >= ads_min."""
    fmt = out_flag(out)
    expr = (f"F_MISSING < {f_missing_max} & "
            f"COUNT(FMT/ADS>={ads_min})/N_SAMPLES >= {sample_frac_min}")
    sh(f"bcftools annotate -x INFO/F_MISSING {q(inp)} -Ou "
       f"| bcftools +fill-tags -Ou -- -t AC,AN,AF,F_MISSING "
       f"| bcftools view -i {q(expr)} -O{fmt} -o {q(out)}", tools=("bcftools",))


def maf_filter(inp: str, out: str, *, maf_min: float = 0.01, maf_max: float = 0.99) -> None:
    """Drop rare and near-fixed alleles by a minor-allele-frequency window."""
    fmt = out_flag(out)
    sh(f"bcftools +fill-tags {q(inp)} -Ou -- -t AC,AN,AF "
       f"| bcftools view -q {maf_min} -Q {maf_max} -O{fmt} -o {q(out)}", tools=("bcftools",))


# --- small bcftools query helpers -------------------------------------------

def _count(path: str) -> int:
    p = subprocess.run(f"bcftools view -H {q(path)} | wc -l", shell=True,
                       executable="/bin/bash", stdout=subprocess.PIPE, text=True)
    return int(p.stdout.strip() or 0)


def _samples(path: str) -> list[str]:
    p = subprocess.run(f"bcftools query -l {q(path)}", shell=True,
                       executable="/bin/bash", stdout=subprocess.PIPE, text=True)
    return [s for s in p.stdout.splitlines() if s]


def _per_sample_covered(path: str, ads_min: int) -> dict[str, int]:
    """Per-sample count of loci with FMT/ADS >= ads_min."""
    cmd = (f"bcftools query -f '[%SAMPLE\\t%ADS\\n]' -i 'FMT/ADS>={ads_min}' {q(path)} "
           "| cut -f1 | sort | uniq -c")
    p = subprocess.run(cmd, shell=True, executable="/bin/bash",
                       stdout=subprocess.PIPE, text=True)
    counts: dict[str, int] = {}
    for line in p.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        n, name = line.split(None, 1)
        counts[name] = int(n)
    return counts


def _write_tmp_list(names: list[str]) -> str:
    import tempfile
    fh = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    fh.write("\n".join(names) + "\n")
    fh.close()
    return fh.name
