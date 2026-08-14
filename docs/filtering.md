# VCF filtering & harmonization

Each step is a parameterized command backed by `bcftools`/`bedtools`, reporting
before/after variant counts.

```bash
plasgenomicsutils hard_qc_filter --input in.bcf --output 01.bcf            # QD/MQ/SOR/RankSums, keep PASS
plasgenomicsutils singleton_filter_add_ads --input 01.bcf --output 02.bcf  # drop singletons, add FORMAT/ADS
plasgenomicsutils tandem_repeat_mask --input 02.bcf --output 03.bcf        # --bed defaults to builtin:pf3d7_tandem_repeats
plasgenomicsutils core_region_filter  --input 03.bcf --output 04.bcf       # keep core genome
plasgenomicsutils paralog_mask        --input 04.bcf --output 05.bcf       # drop paralog/multigene families (optional, see Pipeline)
plasgenomicsutils filter_ad_regenotype --input-vcf 05.bcf --output-vcf 06.bcf  # clean low AD, re-genotype
plasgenomicsutils biallelic_snp_filter --input 06.bcf --output 07.bcf      # keep biallelic SNPs
plasgenomicsutils sample_coverage_filter --input 07.bcf --output 08.bcf
plasgenomicsutils locus_missingness_filter --input 08.bcf --output 09.bcf
plasgenomicsutils maf_filter --input 09.bcf --output 10.bcf --maf-min 0.02 --maf-max 0.98
```

Region masks (`tandem_repeat_mask`, `core_region_filter`, `paralog_mask`) take `--bed`,
a plain path or a bundled asset via `builtin:<name>` (`pf3d7_core_regions`,
`pf3d7_paralog_genes`, `pf3d7_tandem_repeats`).

### Whitelisting regions from a region filter

All three take `--keep-bed`: a BED of regions to keep whatever that filter says. Use it when a
mask is right in general but wrong somewhere specific — a few positions inside a tandem repeat
that are known to be real, or a gene the paralog list catches but you trust.

```bash
plasgenomicsutils tandem_repeat_mask --input 02.bcf --output 03.bcf \
  --keep-bed keep_these.bed
```

It is region-based, so one line covers a whole gene without enumerating its variants:

```
Pf3D7_05_v3	958000	962200	pfmdr1        # every variant in the gene survives this filter
Pf3D7_07_v3	403623	403626	pfcrt-AA76    # or just one codon
```

Three things worth knowing:

- **This exempts a variant from that one rule only.** Whitelisted variants still face every
  other filter in the chain — QC, re-genotyping, coverage, MAF. It is not a "keep no matter
  what" list.
- **The BED is 0-based half-open**, like any BED, and unlike a VCF `POS`. Getting that wrong
  shifts everything a base, so each filter reports how many variants the whitelist actually
  rescued and warns when the answer is zero — a silent no-op is otherwise indistinguishable
  from success.
- **A whitelisted region rescues whole records.** Overlapping any part of a record is enough,
  so a one-base entry keeps a variant whose `REF` is longer than one base (`bedtools` sizes a
  record as `[POS-1, POS-1+len(REF))`, which reaches past a single position).

In the pipeline config it is just another param:

```json
{ "name": "tandem_repeat_mask",
  "params": { "bed": "builtin:pf3d7_tandem_repeats", "keep_bed": "keep_these.bed" } }
```

#### Building one from amino-acid positions

Residue numbers are usually how these exceptions are known ("*pfpx1* 1701 and 1705 are fine").
`aa_intervals()` in `plasgenomicsutilsR` converts them, and its default output is already
0-based half-open — the BED convention — so it can be written straight out:

```r
cds <- read_gff_cds("PlasmoDB-68_Pfalciparum3D7.gff")
codons <- aa_intervals(data.frame(transcript_id = "pfpx1", aa_position = c(1701, 1705)), cds)
write.table(codons[, c("chrom", "start", "end", "name")], "keep_these.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

Do not pass `one_based_output = TRUE` here — that is for numbers you quote in text, and would
put every interval one base off. If a codon straddles an intron its `start`/`end` span the
intron too (`spans_intron` flags it), which for a whitelist only means exempting a little more
than the three bases; use `codon_positions` if you want exactly the coding bases.

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

Each step writes `filtered/NN_<name>.bcf` (indexed) plus a `variant_counts.tsv` tally
(`step`, `kind`, `count`, `path` — `kind` separates a filter's variant count from a report's
row count and from a step that was switched off).

### Turning steps on and off

A step with `"enabled": false` stays in the config and out of the run. JSON has no comments,
so that is how an optional step stays discoverable instead of being something you have to
know exists:

```json
{ "name": "paralog_mask", "enabled": false,
  "params": { "bed": "builtin:pf3d7_paralog_genes" } }
```

**`paralog_mask` ships off.** `core_region_filter` has already removed the subtelomeric
multigene families that mismap worst; most of what `paralog_mask` would take next sits in the
core and is single-copy. Some of that genuinely misbehaves, but a lot of it is fine, and
dropping all of it costs real signal — paralogy is not by itself evidence a locus is
unusable. Set `"enabled": true` when mismapping is specifically what you are controlling for,
or run the step standalone on the region you doubt.

### Why a sample was dropped

`sample_coverage_filter` writes the coverage table its decision is read from, beside the step
that produced it — `filtered/09_sample_coverage_filter_cov_info.tsv`:

| sample | n_loci | n_covered | frac_covered | n_missing_ads | mean_ads | ads_min | frac_min | margin | dropped |
|---|---|---|---|---|---|---|---|---|---|
| RCN13010 | 1901 | 843 | 0.443 | 0 | 15.64 | 5 | 0.5 | −0.057 | True |
| RCN13048 | 1901 | 912 | 0.480 | 0 | 57.26 | 5 | 0.5 | −0.020 | True |
| RCN13065 | 1901 | 964 | 0.507 | 0 | 32.15 | 5 | 0.5 | +0.007 | False |

Sorted worst-covered first, so a surprising drop is the first thing you read. `margin` is
`frac_covered - frac_min`: negative means dropped, and near zero means the sample sat on the
threshold rather than being obviously bad — `RCN13048` above missed by 0.02 while carrying a
mean ADS of 57, i.e. deeply sequenced where it *is* covered and simply absent elsewhere.
That is the distinction `mean_ads` and `n_missing_ads` exist to make, and it usually decides
whether to lower `frac_min` or drop the sample. The same rows go to the console for anything
dropped or within 0.05 of the threshold, and `--cov-table` writes them from the standalone
command.

The table is derived from the same counts as the keep/drop decision, so it always accounts for
what happened rather than being a second measurement that might disagree.

The final callset's **SNP panel BED** is written automatically next to the last step
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

It reports what it did per file — records read, ALT alleles zeroed and removed, records
reduced to ref-only, then records written, records that gained ALTs, and **union sites absent
from the file**.

### Sites one cohort never called

That last count is the one to read. Harmonizing settles which *alleles* the files agree on, not
which *sites* each one contains: a site only cohort A called is simply not written for cohort B.
`bcftools merge` then fills B's samples there with a missing genotype **and a missing
`FORMAT/AD`**, whole-cohort at a time.

That breaks any tool reading `AD` as an integer. hmmibd-rs stops at the first one with:

```
Error: error in processing bcf file: BcfReaderError(NumericaValueEmptyInt)
```

To keep only the sites every cohort called:

```bash
bcftools view -e 'FMT/AD="."' merged.bcf -Ob -o merged.nomiss.bcf
```

Use `FMT/AD="."`, not `FMT/AD[0]="."` — indexing a `Number=R` field per sample does not test
for a missing value and silently matches nothing. Note this is not the same as filtering on
missing genotypes: plenty of sites carry genuine per-sample missing `GT` while `AD` is present,
so `-i 'N_PASS(GT="mis")=0'` would throw away far more than necessary.

The alternative, if you would rather keep those sites, is a caller that tolerates missing depth
— `locus_missingness_filter` will trim them by missingness fraction, but any site absent from a
whole cohort fails that too.

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
