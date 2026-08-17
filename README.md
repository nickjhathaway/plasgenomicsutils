# plasgenomicsutils

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![tests](https://github.com/nickjhathaway/plasgenomicsutils/actions/workflows/tests.yml/badge.svg)](https://github.com/nickjhathaway/plasgenomicsutils/actions/workflows/tests.yml)
<!-- badges: end -->

> **Version 0.2.3** — early development; APIs, defaults, and outputs may change
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

# 2. per-pair / per-SNP / per-group / per-chromosome summaries
#    (--group-col names the metadata column to group samples by; it need not be geographic)
plasgenomicsutils analyze_ibd_matrix --matrix ibd_matrix --meta meta.csv \
    --group-col region --pairwise-group-snp --output ibd_analysis

# 3. global + per-group allele frequencies (single pass over the BCF)
plasgenomicsutils compute_allele_freqs --bcf clean.bcf --meta meta.tsv \
    --group-col region --output afs/

# 4. IBD-based selection statistic (XiR,s), genome-wide and per-group
#    reports Bonferroni AND Benjamini-Hochberg FDR, plus the genomic inflation factor;
#    add --permute 200 when that factor is far from 1: it replaces the assumed chi2(1)
#    with a null drawn from the data, giving a family-wise threshold plus calibrated
#    per-SNP p-values and an FDR over them
plasgenomicsutils ibd_selection_statistic --matrix ibd_matrix \
    --af afs/allele_freqs.tsv.gz --af-group afs/group_allele_freqs.tsv.gz \
    --meta meta.csv --group-col region --output ibd_selection

# per-pair IBD fraction (callable-genome denominator) + SNP density
plasgenomicsutils ibd_fraction_and_snp_density --blocks blocks.hmm.txt \
    --snps snps.bed --snp-format bed --reference pf3d7 --output ibd_frac

# 5. per-gene IBD-block overlap between groups (fraction of pairs whose IBD block
#    overlaps each gene) -> feeds the R gene triangles
plasgenomicsutils ibd_gene_overlap --blocks blocks.hmm.txt --genes genes.tsv \
    --meta meta.csv --group-col region --output gene_overlap.tsv.gz

# 6. which sample PAIRS are IBD over each gene, how much of it they share, and the
#    single-linkage cluster each pair belongs to at that gene (gene_cluster_id)
plasgenomicsutils ibd_gene_pairs --blocks blocks.hmm.txt --genes genes.tsv \
    --output gene_pairs.tsv.gz
```

**Short IBD segments are dropped by default.** Every tool that reads hmm blocks discards
segments with fewer than 15 SNPs or shorter than 15 kb (`--min-block-snp` /
`--min-block-kb`, `0` to disable), matching R's
`ibd_results(min_block_snp = 15, min_block_kb = 15)`. Small blocks are commonly spurious,
and the filter applies to the IBD *evidence* only — the set of pairs that were compared,
the denominator in every fraction, still comes from every row of the blocks file.

**Coordinates are 0-based throughout.** Intervals are half-open `[start, end)` (BED), and
SNP ids are `chr:pos0`, built in exactly one place from `(chrom, pos0)`. An id already
present in an input — a BED name column, or a VCF `ID` set by `bcftools annotate --set-id`
— is never adopted as a key, since whether it used `%POS` or `%POS0` is unknowable from the
file; it is carried alongside as `source_id`. VCF `POS` and `hmmibd-rs` block ends are
converted at the boundary and their numbering never propagates inward. Pass
`--with-pos-vcf` to `compute_allele_freqs` if you also want the 1-based position for
looking variants up by eye.

Outputs carrying SNP labels record the convention (`#snp_coord_system=0-based`) and the
readers verify it, so a table written by an older version is rejected rather than silently
mixed — regenerate the matrix and allele frequencies together.

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

**Stale genotype-linked fields.** If an upstream caller wrote a genotype ploidy that disagrees with a
`Number=G` FORMAT field — e.g. a diploid `GT` forced over hexaploid calls leaves a `PL` of the wrong
length — `bcftools view --trim-alt-alleles` aborts ("Unexpected number of values in FORMAT/PL …").
`strip_stale_format` fixes it: by default (`--mode mismatch`) it nulls such a field **only on the
records where its length is inconsistent** with the genotypes, keeping valid likelihoods elsewhere
(`--mode always` drops the field entirely; default field is `PL`, override with `--fields`).
`filter_ad_regenotype` also does this automatically — it re-genotypes to a fresh `GT`, so it blanks
the now-stale `PL`/`GL` to a consistent length as it writes — and `biallelic_snp_filter` applies the
same surgical fix before trimming, so the default chain never chokes on these fields.

`filter_ad_regenotype` keeps the conventional diploid coding (`0/1` = mixed infection) by default;
`--ploidy {1,2}` sets it explicitly and is validated against the input ploidy per record (greater than
the input errors — genotypes can't be promoted; less warns and trims). Use `--ploidy 1` for haploid
calls.

```bash
plasgenomicsutils strip_stale_format --input calls.bcf --output clean.bcf          # null inconsistent PL
plasgenomicsutils strip_stale_format --input calls.bcf --output clean.bcf --mode always --fields PL GL
```

Or run the whole chain from a JSON config, with a per-step count tally:

```bash
plasgenomicsutils filter_pipeline --emit-default-config pipeline.json   # write a template
plasgenomicsutils filter_pipeline --input in.bcf --config pipeline.json --outdir filtered/
```

Steps write indexed BCFs plus a `variant_counts.tsv`, and the final callset's **SNP-panel
BED** (`filtered/NN_<last>.snps.bed`) is written automatically — it drops straight into
`build_ibd_matrix --snps … --snp-format bed` for the IBD analysis (`--no-snp-bed` to skip).

`maf_filter` can also filter **per group**: pass `--meta samples.tsv --group-col country`
(or set them in the step's `params`) to keep a site if its minor-allele frequency is ≥
`--maf-min` in *any* group. It picks the sites on the combined VCF and applies the union
back to the original, so a variant polymorphic in one country but rare/absent in another is
kept with **all genotypes preserved** — a carrier below the threshold within its own country
keeps its real `0/1`/`1/1` call, with no split-and-merge that would blank it. (A variant
rare in *every* group is dropped.)

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

## Auto completion

To enable bash tab-completion for `plasgenomicsutils` (command names and each command's
options), append the generated script to your `~/.bash_completion` and source it:

```bash
plasgenomicsutils --bash-completion >> ~/.bash_completion
source ~/.bash_completion
```

The same script is also committed at [`etc/bash_completion`](etc/bash_completion) if you
prefer to source a file directly. Completions are queried live from the installed CLI, so
they stay in sync as commands and options change.

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

## Coverage QC (from BAMs)

`coverage_depth_stats` summarises sequencing depth per sample straight from indexed
BAM/CRAM — mean, median, SD, quartiles and the *breadth* at a set of thresholds, per
chromosome and genome-wide. Restrict it to a BED (the core genome is bundled) so
subtelomeric and hypervariable regions do not drag every statistic down:

```bash
plasgenomicsutils coverage_depth_stats \
  --bam-list bams.txt --regions builtin:pf3d7_core_regions \
  --thresholds 1,5,10,20 --window 1000 \
  --output coverage_by_sample.tsv.gz --windows-output coverage_windows.tsv.gz
```

Depth comes from **pysam** by default and from **mosdepth** when it is on `PATH` (it is in
`environment.yml`), which is far faster on whole-genome BAMs — roughly a minute a sample
against four. **The two do not define depth identically**, so the engine used is written
into every output row and should be held fixed across a cohort:

| engine | counts | agrees with |
| --- | --- | --- |
| `pysam` | reads — overlapping mates counted twice | `samtools depth`, base for base |
| `mosdepth` | **fragments** — an overlapping mate pair counts once | runs 2–3% below `pysam` on real Pf WGS |

Fragment depth is the better measure of independent evidence — two mates of one molecule
are one observation, not two — which is the usual reason to prefer mosdepth over
`samtools depth` in the first place. Read depth is what most tools report. Neither is
wrong; mixing them across a cohort is, so the engine used is written into every row.

```bash
# pin fragment depth explicitly rather than depending on what is installed
plasgenomicsutils coverage_depth_stats --bam-list bams.txt --engine mosdepth \
  --regions builtin:pf3d7_core_regions --jobs 8 --output coverage.tsv.gz

# read-level depth, comparable with samtools depth (and the only engine that can
# apply a base-quality floor)
plasgenomicsutils coverage_depth_stats --bam-list bams.txt --engine pysam \
  --min-baseq 13 --jobs 8 --output coverage.reads.tsv.gz
```

`--jobs N` spreads samples across cores. Name your samples in the second column of
`--bam-list` when the BAM filename is not the sample id:

```
/path/4089106922.sorted.bam	IMH07_4089106922
```

`coverage_dropout_regions` then answers the cross-sample question. Selective whole-genome
amplification does not amplify uniformly, and a region that no sample amplifies reads as
*invariant* rather than as missing. This finds the windows below depth in nearly every
sample and merges them into regions:

```bash
plasgenomicsutils coverage_dropout_regions \
  --windows coverage_windows.tsv.gz --regions builtin:pf3d7_core_regions \
  --min-depth 5 --min-frac-samples 0.9 --merge-gap 1000 \
  --genes genes.tsv --output dropouts.tsv.gz --bed-output dropouts.bed
```

Restricting to the core matters here: subtelomeric dropout is expected and would bury the
regions worth acting on. The BED can be fed straight back in as a mask.

Plot both with the R package — `coverage_qc()`, `plot_coverage_summary()`,
`plot_coverage_by_chrom()`, `plot_coverage_dropout()`.

## Sample QC: singleton counts

`singleton_counts` counts, per sample, the variants where it is the only non-reference
carrier. A sample far above the cohort's rate is usually contaminated, mixed-species or
mis-aligned rather than interesting — MalariaGEN drop samples on this criterion when
assembling the Pf analysis sets. Outliers are called by median absolute deviation, so the
outliers being looked for cannot inflate the spread they are measured against:

```bash
plasgenomicsutils singleton_counts --vcf cohort.snps.bcf \
  --min-depth 5 --max-missing-frac 0.2 --mad-cutoff 5 --output singletons.tsv.gz
```

The flag is two-sided, because the tails mean different things:

- **Excess** private variants — contamination, a mixed-species infection, mis-alignment.
- **Deficit** — another sample is absorbing them. The same scan records which partner each
  sample shares its *doubletons* with, and names pairs above `--duplicate-frac` (0.9) as
  near-identical. It deliberately does not call them duplicates: the same parasite
  sequenced twice and one clone infecting two hosts look identical here, and in a
  low-transmission setting the second is the common answer. Check IBD — a clone pair sits
  near IBD 1 — and the collection records before dropping anything.

`--min-depth` (default 5) matters on gVCF-derived callsets, which emit `0/0` at sites with
**no reads at all** rather than `./.`; without it those count as confident reference calls.
Singleton status is judged **within the samples analysed**, so `--samples` changes what
counts as private.

This runs inside `filter_pipeline` too, as a *report* step — it writes a table and passes
the callset through untouched. Position is not arbitrary: `singleton_filter_add_ads` drops
exactly the variants being counted, so the default config puts the report before it (and
warns if a custom config puts it after, where every sample scores zero):

```json
{"steps": [
  {"name": "hard_qc_filter"},
  {"name": "singleton_counts", "report": true, "ext": "tsv",
   "params": {"min_depth": 5, "max_missing_frac": 0.2}},
  {"name": "singleton_filter_add_ads"}
]}
```

## Development

```bash
git clone https://github.com/nickjhathaway/plasgenomicsutils
cd plasgenomicsutils
pip install -e ".[dev]"
pytest
```

## License

GPL-3.
