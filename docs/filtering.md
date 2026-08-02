# VCF filtering & harmonization

Each step is a parameterized command backed by `bcftools`/`bedtools`, reporting
before/after variant counts.

```bash
plasgenomicsutils hard_qc_filter --input in.bcf --output 01.bcf            # QD/MQ/SOR/RankSums, keep PASS
plasgenomicsutils singleton_filter_add_ads --input 01.bcf --output 02.bcf  # drop singletons, add FORMAT/ADS
plasgenomicsutils tandem_repeat_mask --input 02.bcf --output 03.bcf        # --bed defaults to builtin:pf3d7_tandem_repeats
plasgenomicsutils core_region_filter  --input 03.bcf --output 04.bcf       # keep core genome
plasgenomicsutils paralog_mask        --input 04.bcf --output 05.bcf       # drop paralog/multigene families
plasgenomicsutils filter_ad_regenotype --input-vcf 05.bcf --output-vcf 06.bcf  # clean low AD, re-genotype
plasgenomicsutils biallelic_snp_filter --input 06.bcf --output 07.bcf      # keep biallelic SNPs
plasgenomicsutils sample_coverage_filter --input 07.bcf --output 08.bcf
plasgenomicsutils locus_missingness_filter --input 08.bcf --output 09.bcf
plasgenomicsutils maf_filter --input 09.bcf --output 10.bcf --maf-min 0.02 --maf-max 0.98
```

Region masks (`tandem_repeat_mask`, `core_region_filter`, `paralog_mask`) take `--bed`,
a plain path or a bundled asset via `builtin:<name>` (`pf3d7_core_regions`,
`pf3d7_paralog_genes`, `pf3d7_tandem_repeats`).

!!! note "Step order depends on your input"
    For an already-biallelic core-SNP callset the default order works. For a raw,
    multiallelic joint callset, run `biallelic_snp_filter` **before**
    `singleton_filter_add_ads` (whose `FORMAT/ADS` and the downstream MAF assume
    biallelic sites). MAF filtering is per-population — run once per cohort.

## AD cleaning and re-genotyping

`filter_ad_regenotype` judges depth by `AD`/`ADS` (not `DP`), zeros sub-threshold
allele depths per sample, recomputes `ADS`, and re-genotypes conservatively: alleles
are ranked by depth and a heterozygote is called only when the minor allele's
within-sample frequency is ≥ `--het-min-af` (default 0.2). This preserves the
mixed-infection signal Fws/COI analyses depend on. `--restrict-to-called-alleles`
narrows the caller's existing genotypes instead of re-deriving from AD.

**Ploidy.** The default keeps the conventional diploid coding used for *Plasmodium*
(`0/1` = mixed infection). `--ploidy {1,2}` sets it explicitly and is validated against
the input ploidy per record: a value **greater** than the input errors (genotypes and
likelihoods cannot be promoted), **equal** is fine, and **less** warns and trims. Use
`--ploidy 1` for true haploid calls (the single best-supported allele).

**Stale genotype-linked fields.** Re-genotyping writes a fresh `GT`, so a caller's
`Number=G` fields (`PL`/`GL`) from a different ploidy — e.g. a diploid `GT` forced over
hexaploid calls — become inconsistent and would break
`bcftools view --trim-alt-alleles`. `filter_ad_regenotype` blanks them to a consistent
length automatically as it writes. For files that did **not** go through re-genotyping,
`strip_stale_format` does the same on its own:

```bash
plasgenomicsutils strip_stale_format --input calls.bcf --output clean.bcf            # null only the inconsistent PL records (surgical, keeps valid ones)
plasgenomicsutils strip_stale_format --input calls.bcf --output clean.bcf --mode always --fields PL GL   # drop the fields entirely
```

`biallelic_snp_filter` applies the same surgical fix before trimming, so the default
chain never chokes on these fields.

`maf_filter` takes `--maf-min`/`--maf-max`; since the bounds are usually symmetric,
`--maf-max` defaults to `1 - maf_min` when unset (so `maf_min = 0.02` gives a `[0.02, 0.98]`
window). Set both for an asymmetric window.

**Per-group MAF.** With `--meta <table.tsv> --group-col country` (a per-sample metadata
table with a `sample` column and the group column), a site is kept if its minor-allele
frequency is ≥ `--maf-min` in **any** group. It is computed on the combined VCF — the
per-group frequencies pick the *sites* to keep, then that union is applied back to the
original, so **every sample's genotypes are preserved**. A carrier whose variant is below
the threshold *within its own group* keeps its real `0/1`/`1/1` call whenever the site is
kept via another group (and a `0/0` sample stays `0/0`) — there is no split-and-merge, so
nothing is blanked. A variant that is rare in *every* group is dropped.

```bash
plasgenomicsutils maf_filter --input in.bcf --output maf.bcf \
  --maf-min 0.02 --meta samples.tsv --group-col country
```

## Pipeline

Run an ordered, config-driven chain and tally counts per step:

```bash
plasgenomicsutils filter_pipeline --emit-default-config pipeline.json
plasgenomicsutils filter_pipeline --input in.bcf --config pipeline.json --outdir filtered/
```

Each step writes `filtered/NN_<name>.bcf` (indexed) plus a `variant_counts.tsv` tally. The
final callset's **SNP panel BED** is written automatically next to the last step
(`filtered/NN_<last>.snps.bed`) — this is the panel the IBD tools read, so it feeds straight
into `build_ibd_matrix`:

```bash
plasgenomicsutils build_ibd_matrix --blocks hmmibd.tsv.gz \
  --snps filtered/10_maf_filter.snps.bed --snp-format bed --output ibd_matrix.npz
```

The BED's `chrom:pos` name column matches the SNP ids you'd get from the VCF, so a panel
taken as BED or VCF is interchangeable. Pass `--no-snp-bed` to skip it.

## Harmonization

Align the allele sets of separately-called cohorts before `bcftools merge`:

```bash
plasgenomicsutils harmonize_bcf --files cohortA.bcf cohortB.bcf cohortC.bcf --stub harmonized
```

It strips stale `Number=A/R/G` fields that a reshaped `AD` would otherwise break, builds
the per-site ALT union, and re-genotypes only changed samples. Indel-context records are
dropped by default (`--keep-indels` to disable).

## Strand-bias artifacts

Illumina sequence-specific errors produce fake heterozygous calls whose ALT reads sit
almost entirely on one strand. Scan a biallelic callset carrying per-strand depths and
emit a blacklist BED (feed it to `tandem_repeat_mask --bed`):

```bash
bcftools mpileup -a FORMAT/ADF,FORMAT/ADR -f ref.fa -b bams.txt -Ou \
  | bcftools norm -m- -f ref.fa -Ob -o strand.bcf
plasgenomicsutils strand_bias_scan --input-vcf strand.bcf --out-tsv verdicts.tsv --out-bed sse.bed
```

Characterize one site at read level (and optionally dump ALT reads):

```bash
plasgenomicsutils strand_read_check --bam sample.bam --pos Pf3D7_12_v3:975431 \
  --ref ref.fa --alt-base A --extract-reads
```

See [Strand-bias artifacts](strand_bias_artifact_exclusion.md) for the detection rules.
