# Coverage QC

Two commands read the BAMs and write tables; `plasgenomicsutilsR` reads and plots them
(see its *Coverage QC* article).

```bash
# per-sample, per-chromosome depth, plus genome-wide rows
plasgenomicsutils coverage_depth_stats --bam-list bams.txt --thresholds 10 \
  --output coverage_by_sample.tsv.gz --windows-output coverage_windows.tsv.gz

# windows below depth in nearly every sample, merged into regions
plasgenomicsutils coverage_dropout_regions --windows coverage_windows.tsv.gz \
  --min-depth 5 --min-frac-samples 0.9 --output dropout_regions.tsv \
  --bed-output dropout_regions.bed
```

## Breadth matters more than mean depth

Selective whole-genome amplification can give a respectable average while leaving much of
the genome at zero, and only the breadth column shows that. `--thresholds 10` writes
`pct_ge_10x`, which is the column to read next to `mean`. On one real cohort the sample
that failed QC had a mean of 122x — and a median of 33x, with a quarter of the core genome
under 10x.

## The engines do not agree, by definition

`--engine mosdepth` counts **fragments** (an overlapping mate pair once); `--engine pysam`
counts **reads**. That puts mosdepth a couple of percent lower everywhere. Fragment depth
is the better measure of independent evidence; either is fine as long as a cohort uses one
throughout, and the `engine` column records which was used. `--engine auto` (default)
takes mosdepth when it is on the PATH.

## Dropouts are a cross-sample question

A window empty in *one* sample is missing data for that sample. A window empty in
*everyone* is not missing data at all — it silently reads as invariant, and every
downstream statistic treats it as a region where nobody carries an alternate allele.
`coverage_dropout_regions` is what finds those, and `--bed-output` writes them in a form
the region filters accept:

```bash
plasgenomicsutils core_region_filter --input clean.bcf --output cleaner.bcf \
  --bed builtin:pf3d7_core_regions --keep-bed resistance_loci.bed
```
