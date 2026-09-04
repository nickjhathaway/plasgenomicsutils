# Command reference

Everything runs through one runner; commands share a flat, globally unique namespace.

```bash
plasgenomicsutils --list           # this catalog
plasgenomicsutils <command> -h     # options for one command
```

## `fws`

| Command | Description |
|---|---|
| `calculate_fws` | Per-sample Fws within-host diversity (`moimix::getFws`) from a VCF or AD table |

## `ibd`

| Command | Description |
|---|---|
| `build_ibd_matrix` | Build binary (pairs × SNPs) IBD matrix from hmmibd-rs blocks |
| `compute_allele_freqs` | Compute global + per-region allele frequencies (single pass) |
| `analyze_ibd_matrix` | Per-pair / per-SNP / per-region / per-chromosome IBD summaries |
| `ibd_selection_statistic` | IBD-based selection statistic (XiR,s), genome-wide and per-region |
| `ibd_fraction_and_snp_density` | Per-pair IBD fraction (callable denominator) and SNP density |
| `ibd_gene_overlap` | Fraction of pairs whose IBD block overlaps each gene, per group pair |
| `ibd_gene_pairs` | Sample pairs IBD over each gene, with how much of the gene is covered |

## `ld`

| command | what it does |
| --- | --- |
| `ld_decay` | Mean r² vs SNP-pair distance per group: how fast LD decays |

## `coverage`

| command | what it does |
| --- | --- |
| `coverage_depth_stats` | Per-sample depth: mean/median/SD and breadth at thresholds, per chromosome |
| `coverage_dropout_regions` | Regions below depth in nearly every sample (sWGA amplification dropouts) |

## `vcf`

Making a callset, merging one, converting one.

| Command | Description |
|---|---|
| `call_variants` | Call variants with bcftools, annotated for `hard_qc_filter`, parallel over regions |
| `harmonize_bcf` | Harmonize ALT sets of separately-called cohorts for `bcftools merge` |
| `vcf_to_bed` | Convert a VCF/BCF to 0-based BED (stdout by default) |

## `vcf_filter_pipeline`

`filter_pipeline` runs these as an ordered chain; each is also a standalone command. They are listed in the order the default config runs them rather than alphabetically — the order is the point, and a step reads differently depending on what has already happened to the callset.

| Command | Description |
|---|---|
| `filter_pipeline` | Run an ordered, config-driven chain of filtering steps, tallying counts |
| `no_alt_filter` | Drop records with no ALT allele (non-variant positions), counted separately |
| `hard_qc_filter` | GATK-style hard filter on INFO metrics (QD/MQ/SOR/RankSums), keep PASS |
| `singleton_filter_add_ads` | Drop near-private variants and add the `FORMAT/ADS` summed-depth tag |
| `tandem_repeat_mask` | Remove variants overlapping a tandem-repeat BED |
| `core_region_filter` | Keep only variants inside the core-genome BED |
| `paralog_mask` | Remove variants overlapping paralogous/multigene-family genes |
| `filter_ad_regenotype` | Clean within-sample AD artifacts by depth/frequency, then re-genotype |
| `biallelic_snp_filter` | Keep biallelic SNPs, trimming ALT alleles unused after re-genotyping |
| `sample_coverage_filter` | Drop low-coverage samples; refresh AC/AN/AF |
| `locus_missingness_filter` | Keep loci with low missingness and high per-sample coverage |
| `maf_filter` | Keep variants within a minor-allele-frequency window |
| `fws_filter` | Keep only monoclonal samples (Fws >= a threshold); refresh AC/AN/AF |
| `strip_stale_format` | Strip stale genotype-linked FORMAT fields (e.g. `PL`) that no longer match the genotypes |

## `vcf_reporting`

Reads a callset and writes a table; never changes the data.

| Command | Description |
|---|---|
| `singleton_counts` | Per-sample count of variants where it is the only non-reference carrier |
| `strand_bias_scan` | Flag strand-bias (SSE) fake-het artifacts from `FORMAT/ADF+ADR`; emit a blacklist BED |
| `strand_read_check` | Read-level strand-bias diagnostic at one site (+ optional ALT-read extraction) |
| `variant_spacing` | Per-chromosome gaps between consecutive variants, and density per cM |
| `wsaf_profile` | Per-sample: is there a dominant clone, and what `--min-freq` reduces the sample to it |
