# plasgenomicsutils

A collection of utilities for **post processing Plasmodium genomics data** —
VCF/BCF filtering and harmonization, IBD (identity-by-descent) post-analysis, and
(planned) Fws within-host diversity.

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

## Development

```bash
git clone https://github.com/nickjhathaway/plasgenomicsutils
cd plasgenomicsutils
pip install -e ".[dev]"
pytest
```

## License

GPL-3.
