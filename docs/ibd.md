# IBD post-analysis

Downstream of the `hmmibd-rs` IBD caller, these commands turn IBD blocks into a binary
matrix and per-pair / per-SNP / per-region summaries, allele frequencies, and an
EIGENSTRAT-style selection statistic. Their tables feed the R plots in
[plasgenomicsutilsR](https://nickjhathaway.github.io/plasgenomicsutilsR/).

```bash
# hmmibd-rs blocks + a SNP panel -> sparse (pairs x SNPs) binary matrix
plasgenomicsutils build_ibd_matrix --blocks blocks.tsv --snps panel.vcf.gz --out ibd_matrix

# per-pair / per-SNP / per-region / per-chromosome summaries
plasgenomicsutils analyze_ibd_matrix --matrix ibd_matrix.npz --metadata samples.tsv \
  --region-col region --out ibd_analysis

# population + per-region allele frequencies (single pass over the BCF)
plasgenomicsutils compute_allele_freqs --input cohort.bcf --regions samples.tsv --out allele_freqs.tsv.gz

# IBD selection statistic (XiR,s), genome-wide and per region
plasgenomicsutils ibd_selection_statistic --matrix ibd_matrix.npz \
  --allele-freqs allele_freqs.tsv.gz --out ibd_selection_analysis

# per-pair IBD fraction (callable-length denominator) and SNP density
plasgenomicsutils ibd_fraction_and_snp_density --blocks blocks.tsv --snps panel.vcf.gz \
  --reference pf3d7 --out ibd_fraction
```

Block-size thresholds, the region column, and per-analysis parameters are flags rather
than hard-coded, so a single tool serves different cohorts. The SNP panel and BED
derivation use `bcftools query` throughout.

Run `plasgenomicsutils <command> -h` for the full option list of any step.
