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

# which sample PAIRS are IBD over each gene, and how much of the gene each block covers
plasgenomicsutils ibd_gene_pairs --blocks blocks.tsv.gz --genes genes.tsv \
  --output gene_pairs.tsv.gz
```

## Short IBD segments are dropped by default

Small IBD blocks are commonly spurious, so every tool that reads hmm blocks discards
segments with fewer than **15 SNPs** or shorter than **15 kb** — `--min-block-snp` /
`--min-block-kb`, `0` to disable either. The filter is built in rather than something you
have to apply beforehand, and it is applied to the IBD **evidence** only: the set of pairs
that were compared (the denominator in every fraction) still comes from every row of the
blocks file, so a pair whose only segment was short still counts as compared and simply
contributes no sharing.

```bash
# stricter than the default
plasgenomicsutils build_ibd_matrix --blocks blocks.tsv.gz --snps panel.snps.bed \
  --snp-format bed --min-block-snp 20 --min-block-kb 25 --output ibd_matrix

# keep every segment (what you would get by pre-filtering nothing)
plasgenomicsutils build_ibd_matrix --blocks blocks.tsv.gz --snps panel.snps.bed \
  --snp-format bed --min-block-snp 0 --min-block-kb 0 --output ibd_matrix_unfiltered
```

`build_ibd_matrix`, `ibd_gene_overlap`, `ibd_gene_pairs` and
`ibd_fraction_and_snp_density` all take the flags, with the same defaults as R's
`ibd_results(min_block_snp = 15, min_block_kb = 15)`, so the two sides agree. The filter
is not free: on a typical cohort it removes a large share of the raw segments, which shifts
every fraction downward — regenerate the matrix and any per-gene tables together.

`ibd_gene_overlap` reads the hmm blocks directly (a pair counts when any IBD segment
overlaps the gene, even with no panel SNP inside), keeps only IBD segments for the
numerator, and uses **every** compared pair as the denominator (so never-IBD pairs still
count). The `--genes` table needs `name, chr` (or `chrom`), `start`, `end`. Feed the output
to R's `ibd_results(gene_overlap = ...)`, or load the blocks in R
(`ibd_results(blocks = ..., meta = ...)`) to compute the same thing for ad-hoc genes.

## Which pairs share a gene

`ibd_gene_overlap` gives the *fraction* per group pair; `ibd_gene_pairs` gives the
underlying adjacency list — one row per sample pair x IBD block x gene. Pairs with no IBD
over a gene are simply absent, and a pair repeats only when it has several separate
segments spanning the gene.

Worked example over the eight drug-resistance genes, starting from the gene table the R
package ships (`PF_EXAMPLE_DRUG_GENES`: pfcrt, pfdhfr, pfmdr1, pfdhps, pfkelch13, pfaat1,
pfgch1, pfpx1). Write it out once from R:

```r
library(plasgenomicsutilsR)
g <- PF_EXAMPLE_DRUG_GENES[, c("name", "Pf3D7_chrom", "start", "end", "gene_id")]
names(g)[2] <- "chr"
write.table(g, "drug_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

...then list the sharing pairs for all eight at once:

```bash
plasgenomicsutils ibd_gene_pairs \
  --blocks final_set_filtered/10_maf_filter_hmm.hmm.txt \
  --genes drug_genes.tsv \
  --output drug_gene_pairs.tsv.gz

# only the pairs sharing a whole gene
plasgenomicsutils ibd_gene_pairs --blocks final_set_filtered/10_maf_filter_hmm.hmm.txt \
  --genes drug_genes.tsv --complete-only --output drug_gene_pairs.complete.tsv.gz

# one gene, and only substantial sharing
plasgenomicsutils ibd_gene_pairs --blocks final_set_filtered/10_maf_filter_hmm.hmm.txt \
  --genes drug_genes.tsv --gene pfcrt --min-percent-covered 50 \
  --output pfcrt_pairs.tsv.gz
```

Output columns: `sample1`, `sample2` (order-normalised), `chr`, `block_start`/`block_end`
(the IBD segment), `gene`, `gene_id`, `gene_start`/`gene_end`, `coverage`
(`complete`/`partial`), `covered_start`/`covered_end` (the gene's own bounds when
complete), `covered_bp` and `percent_covered`. All coordinates are 0-based half-open.

The R equivalent returns a tibble instead of writing a file, and the two agree row for row:

```r
ibd <- ibd_results(blocks = "final_set_filtered/10_maf_filter_hmm.hmm.txt",
                   genes = PF_EXAMPLE_DRUG_GENES)
pairs <- ibd$gene_ibd_pairs()                       # all eight genes in one tibble
subset(pairs, gene == "pfcrt" & coverage == "complete")
```

(The SNP panel can be a VCF/BCF as above, or the `.snps.bed` that `filter_pipeline` writes —
use `--snp-format bed`.)

Block-size thresholds, the grouping column (`--group-col`), and per-analysis parameters are flags rather
than hard-coded, so a single tool serves different cohorts. The SNP panel and BED
derivation use `bcftools query` throughout.

Run `plasgenomicsutils <command> -h` for the full option list of any step.
