# plasgenomicsutils

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![tests](https://github.com/nickjhathaway/plasgenomicsutils/actions/workflows/tests.yml/badge.svg)](https://github.com/nickjhathaway/plasgenomicsutils/actions/workflows/tests.yml)
<!-- badges: end -->

> **Version 0.2.0** — early development; APIs, defaults, and outputs may change
> between versions.

A collection of utilities for **post processing Plasmodium genomics data** —
VCF/BCF filtering and harmonization, IBD (identity-by-descent) post-analysis, and
Fws within-host diversity.

The heavy compute lives here in Python; **visualization lives in the companion R
package [`plasgenomicsutilsR`](https://github.com/nickjhathaway/plasgenomicsutilsR)**.
These utilities are written to be dataset- and *species*-agnostic: the only
species-specific facts (chromosome lengths, genetic-map rate) live in a small
reference registry (`--reference`, default `pf3d7`), so other species can be added
without touching the algorithms.

## Install

```bash
pip install git+https://github.com/nickjhathaway/plasgenomicsutils
```

Requires Python ≥ 3.10 and, on `PATH`, `bcftools`/`bedtools` for the VCF filtering
commands. Python deps (`numpy pandas scipy pysam cyvcf2`) install automatically.

## Usage

One runner, compartmentalized subcommands (run `plasgenomicsutils --list` for the
full catalog, `plasgenomicsutils <command> -h` for a command's options):

```bash
plasgenomicsutils --list
```

### IBD post-analysis

Downstream of an IBD caller (we use [`hmmibd-rs`](https://github.com/bguo068/hmmibd-rs)):

```bash
# 1. binary (pairs x SNPs) IBD matrix from hmmibd-rs blocks + a SNP panel (VCF or BED)
plasgenomicsutils build_ibd_matrix --blocks blocks.hmm.txt --snps snps.bed \
    --snp-format bed --output ibd_matrix

# 2. per-pair / per-SNP / per-region / per-chromosome summaries
plasgenomicsutils analyze_ibd_matrix --matrix ibd_matrix --meta meta.csv \
    --region-col region --pairwise-region-snp --output ibd_analysis

# 3. global + per-region allele frequencies (single pass over the BCF)
plasgenomicsutils compute_allele_freqs --bcf clean.bcf --meta meta.tsv \
    --region-col region --zero-based --output afs/

# 4. IBD-based selection statistic (XiR,s), genome-wide and per-region
plasgenomicsutils ibd_selection_statistic --matrix ibd_matrix \
    --af afs/allele_freqs.tsv.gz --af-region afs/region_allele_freqs.tsv.gz \
    --meta meta.csv --region-col region --output ibd_selection

# per-pair IBD fraction (callable-genome denominator) + SNP density
plasgenomicsutils ibd_fraction_and_snp_density --blocks blocks.hmm.txt \
    --snps snps.bed --snp-format bed --reference pf3d7 --output ibd_frac
```

SNP-panel label coordinates must be consistent between the matrix and the allele
frequencies: build the matrix from a BED (0-based) and run `compute_allele_freqs`
with `--zero-based`, or build from a VCF (1-based labels) and omit `--zero-based`.

### VCF/BCF filtering & harmonization

Individual, parameterized filtering steps (each reports before/after variant
counts), backed by `bcftools`/`bedtools`:

```bash
plasgenomicsutils hard_qc_filter --input in.bcf --output 01.bcf          # QD/MQ/SOR/RankSums, keep PASS
plasgenomicsutils singleton_filter_add_ads --input 01.bcf --output 02.bcf # drop singletons, add FORMAT/ADS
plasgenomicsutils tandem_repeat_mask --input 02.bcf --output 03.bcf      # --bed defaults to builtin:pf3d7_tandem_repeats
plasgenomicsutils core_region_filter  --input 03.bcf --output 04.bcf      # keep core genome (builtin:pf3d7_core_regions)
plasgenomicsutils paralog_mask        --input 04.bcf --output 05.bcf      # drop paralog/multigene families
plasgenomicsutils filter_ad_regenotype --input-vcf 05.bcf --output-vcf 06.bcf  # clean low AD, re-genotype
plasgenomicsutils biallelic_snp_filter --input 06.bcf --output 07.bcf     # keep biallelic SNPs (trims artifact alleles)
plasgenomicsutils sample_coverage_filter --input 07.bcf --output 08.bcf
plasgenomicsutils locus_missingness_filter --input 08.bcf --output 09.bcf
plasgenomicsutils maf_filter --input 09.bcf --output 10.bcf --maf-min 0.02 --maf-max 0.98
```

Region masks (`tandem_repeat_mask`, `core_region_filter`, `paralog_mask`) take `--bed`, which
accepts a plain path or a bundled Pf3D7 asset via `builtin:<name>` (`pf3d7_core_regions`,
`pf3d7_paralog_genes`, `pf3d7_tandem_repeats`) — the defaults. **`biallelic_snp_filter` runs after
`filter_ad_regenotype` on purpose:** re-genotyping zeroes a single sample's low-level artifact allele,
so a site that merely *looked* multiallelic is trimmed back to a genuine biallelic SNP and kept,
rather than being discarded by a naive `-m2 -M2`.

Or run the whole chain from a JSON config, with a per-step count tally:

```bash
plasgenomicsutils filter_pipeline --emit-default-config pipeline.json   # write a template
plasgenomicsutils filter_pipeline --input in.bcf --config pipeline.json --outdir filtered/
```

**Step order depends on your input.** The chain is config-driven, so you pick the
order; two regimes are common:

- *Already-biallelic input* (e.g. a core-SNP callset): the default order works as
  shipped — `singleton_filter_add_ads` can run early because every site is
  biallelic.
- *Raw, multiallelic input* (e.g. a fresh joint callset): put `biallelic_snp_filter`
  **before** `singleton_filter_add_ads`. The latter's `FORMAT/ADS` is a per-sample
  sum over all AD entries, but the singleton test and downstream MAF assume
  biallelic sites, so reduce to biallelic SNPs first.

MAF filtering is per-population — run the pipeline once per cohort/country, not on
a pooled multi-population VCF. Region masks default to bundled Pf3D7 BEDs
(`builtin:pf3d7_core_regions` / `_paralog_genes` / `_tandem_repeats`); pass a path
to `--bed` (or the config `params.bed`) to use your own.

`filter_ad_regenotype` and `harmonize_bcf` re-genotype haploid calls
conservatively: alleles are ranked by depth and a heterozygote is called only
when the minor allele's within-sample frequency is ≥ `--het-min-af` (default
0.2), otherwise the sample is homozygous for the major allele.

By default the genotype is re-derived from the cleaned AD alone, because a
callset's `GT` can disagree with its own read counts and likelihoods. A case seen
in practice: variants called at a higher ploidy to pick up low-abundance minor
clones in polyclonal infections, then run through a haploid-assuming "diploidize"
script that keeps only each genotype's first allele. Because callers write unphased
alleles sorted ascending, every partially-alternate call collapses to hom-ref —
records end up as `GT=0/0` with `AD=8,200` (8 reference vs 200 alternate reads),
while the untouched `PL` still makes the all-reference state astronomically
unlikely. Re-deriving from AD recovers that signal, and `--het-min-af` marks the
intermediate cases heterozygous — which for a haploid organism is the
mixed-infection signal Fws/COI analyses depend on.

Pass `--restrict-to-called-alleles` to instead only narrow the caller's existing
genotype (a call can lose support and go missing, but never gain a new allele).
That mode reproduces upstream-caller genotypes exactly, which is what you want
when re-running an existing analysis for bit-identical results.

Harmonize separately-called cohorts so their allele sets line up before
`bcftools merge` (inputs must be coordinate-sorted):

```bash
plasgenomicsutils harmonize_bcf --files cohortA.bcf cohortB.bcf cohortC.bcf --stub harmonized
# then: bcftools merge --merge snps harmonized_*.vcf ... | bcftools +fill-tags ...
```

Defaults to VCF output (`--output-format v`); `z` (bgzipped VCF) and `b` (BCF)
are produced by writing VCF and converting with `bcftools`. Harmonize never
writes BCF directly: when a record's allele count is reduced during cleaning,
pysam leaves the binary `Number=R` AD array at the old length, which reads back
fine as text but breaks `bcftools merge` — converting from VCF regenerates AD
consistently.

Indel-context records are **dropped from every input by default** (`--keep-indels`
to disable) — including no-ALT `REF>.` records that carry the `INDEL` flag, which
`bcftools view --exclude-type indels` does *not* remove (a `.` ALT isn't typed as
an indel), plus any multi-base REF/ALT. This means you don't have to pre-filter
indels on each cohort.

When more than one record still shares a position, the record carrying real ALT
alleles is kept — so the true SNP genotypes are preserved rather than being
clobbered by an overlapping no-ALT record. Positions where more than one record
carries real ALTs are flagged; normalize those inputs with `bcftools norm -m -any`
first.

Allele-dependent FORMAT fields that harmonize cannot recompute after reshaping
alleles — any `Number=A/R/G` field other than `AD` (which it does maintain),
e.g. `PL` — are stripped from the output, since a stale `PL` would otherwise make
`bcftools merge` fail with a FORMAT length mismatch. `AD`, `GT` and scalar fields
are preserved.

### Strand-bias artifact diagnostics

Illumina sequence-specific errors (SSE) at GC islands, G-homopolymers, and GGC
motifs produce **fake heterozygous calls** whose ALT reads sit almost entirely on
one strand, at low base quality. In copy-number-from-VAF work they land near the
0.25/0.33 fractions a real 3–4 copy locus shows, so they mimic the signal they
corrupt. A real het is strand-balanced; an artifact is strand-restricted.

Batch-scan a biallelic callset that carries per-strand depths for these artifacts,
writing a per-`(site, sample)` verdict table and a blacklist BED of confirmed
positions (feed that BED to `tandem_repeat_mask --bed` to drop them):

```bash
bcftools mpileup -a FORMAT/ADF,FORMAT/ADR -f ref.fa -b bams.txt -Ou \
  | bcftools norm -m- -f ref.fa -Ob -o strand.bcf                      # ADF/ADR, biallelic
plasgenomicsutils strand_bias_scan --input-vcf strand.bcf \
  --out-tsv verdicts.tsv --out-bed sse_blacklist.bed
```

The verdict drops a call when the alt is strand-restricted (minor/major strand-VAF
ratio `< --ratio-hard`, *and* the minor strand had the depth to detect it, so
genuinely shallow strands are protected) or the Fisher strand-bias Phred exceeds
`--sb-hard`. Because SSE is position-reproducible, a position is blacklisted once
`--min-drop-samples` samples flag it.

To characterize one suspicious site at read level (strand, base quality, true
sequencing cycle, soft-clip/chimera geometry, priming-start spike) and optionally
dump the ALT reads for a sequence viewer:

```bash
plasgenomicsutils strand_read_check --bam sample.bam --pos Pf3D7_12_v3:975431 \
  --ref ref.fa --alt-base A --extract-reads
```

The theory, detection signature, and exclusion rules are in
[docs/strand_bias_artifact_exclusion.md](docs/strand_bias_artifact_exclusion.md).

## Fws (within-host diversity)

`calculate_fws` computes the per-sample Fws statistic (Manske 2012) — a monoclonal
infection scores Fws ~ 1, a polyclonal one lower — from per-sample allele depths,
read from **either a VCF/BCF or a bcftools-query AD table**. It reimplements
`moimix::getFws`:

```bash
plasgenomicsutils calculate_fws --input-vcf cohort.snps.bcf --out fws.tsv
# or from an AD table (CHROM POS REF ALT then one "ref,alt" per sample):
plasgenomicsutils calculate_fws --ad-table ad.tsv --samples samples.txt --out fws.tsv
```

The default `--estimator regression` reproduces `moimix::getFws` to full float
precision (validated against moimix on public pf7 data — see
[docs/fws_moimix_reconciliation.md](docs/fws_moimix_reconciliation.md)): Fws = 1 − β,
β the slope of a through-origin regression of per-sample on population heterozygosity
across 10 MAF bins. `--estimator ratio` is a simpler summed-binned-mean estimator
(`1 − Σ mean(Hw)/Σ mean(Hs)`) — the two are not interchangeable, so don't mix a
threshold tuned on one with the other's values. moimix parity uses the defaults
(`--min-depth 0 --min-alt-samples 0`).

`moimix::getFws` uses every biallelic record's `AD` regardless of allele string, so
the VCF reader does too — pass a SNP-filtered callset (or `--snps-only`) if you want
SNPs only. `--population-name` tags every row for later cross-cohort merging;
`--exclude-call-regions` drops CNV windows whose within-sample heterozygosity would
otherwise depress Fws.

## Development

```bash
git clone https://github.com/nickjhathaway/plasgenomicsutils
cd plasgenomicsutils
pip install -e ".[dev]"
pytest
```

## License

GPL-3.
