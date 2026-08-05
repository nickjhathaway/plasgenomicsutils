# IBD post-analysis

Downstream of the `hmmibd-rs` IBD caller, these commands turn IBD blocks into a binary
matrix and per-pair / per-SNP / per-region summaries, allele frequencies, and an
EIGENSTRAT-style selection statistic. Their tables feed the R plots in
[plasgenomicsutilsR](https://nickjhathaway.github.io/plasgenomicsutilsR/).

```bash
# hmmibd-rs blocks + a SNP panel -> sparse (pairs x SNPs) binary matrix (writes ibd_matrix.npz)
plasgenomicsutils build_ibd_matrix --blocks blocks.tsv.gz --snps panel.vcf.gz --snp-format vcf \
  --output ibd_matrix

# per-pair / per-SNP / per-group / per-chromosome summaries
plasgenomicsutils analyze_ibd_matrix --matrix ibd_matrix --meta samples.tsv \
  --group-col region --output ibd_analysis

# global + per-group allele frequencies, single pass over the BCF
# (writes allele_freqs.tsv.gz and group_allele_freqs.tsv.gz in the output directory)
plasgenomicsutils compute_allele_freqs --bcf cohort.bcf --meta samples.tsv \
  --group-col region --output freqs/

# IBD selection statistic (XiR,s), genome-wide and per group
plasgenomicsutils ibd_selection_statistic --matrix ibd_matrix \
  --af freqs/allele_freqs.tsv.gz --af-group freqs/group_allele_freqs.tsv.gz \
  --meta samples.tsv --group-col region --output ibd_selection

# per-pair IBD fraction (callable-length denominator) and SNP density
plasgenomicsutils ibd_fraction_and_snp_density --blocks blocks.tsv.gz --snps panel.vcf.gz \
  --snp-format vcf --reference pf3d7 --output ibd_fraction

# per-gene IBD-block overlap between groups (for the R gene triangles): the fraction of
# pairs whose IBD block overlaps each gene, NOT pairs that are IBD at a SNP inside it
plasgenomicsutils ibd_gene_overlap --blocks blocks.tsv.gz --genes genes.tsv \
  --meta samples.tsv --group-col region --output gene_overlap.tsv.gz
```

`ibd_gene_overlap` reads the hmm blocks directly (a pair counts when any IBD segment
overlaps the gene, even with no panel SNP inside), keeps only IBD segments for the
numerator, and uses **every** compared pair as the denominator (so never-IBD pairs still
count). The `--genes` table needs `name, chr` (or `chrom`), `start`, `end`. Feed the output
to R's `ibd_results(gene_overlap = ...)`, or load the blocks in R
(`ibd_results(blocks = ..., meta = ...)`) to compute the same thing for ad-hoc genes.

(The SNP panel can be a VCF/BCF as above, or the `.snps.bed` that `filter_pipeline` writes —
use `--snp-format bed`.)

Block-size thresholds, the grouping column (`--group-col`), and per-analysis parameters are flags rather
than hard-coded, so a single tool serves different cohorts. The SNP panel and BED
derivation use `bcftools query` throughout.

Run `plasgenomicsutils <command> -h` for the full option list of any step.
