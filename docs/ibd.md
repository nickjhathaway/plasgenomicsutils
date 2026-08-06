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

## Calling a selection peak significant

`ibd_selection_statistic` reports up to four thresholds, all written to
`<output>.*.threshold.txt` and applied as flag columns of the stats table. The last two
appear only with `--permute`:

| column | threshold | controls | rests on |
|---|---|---|---|
| `significant` | Bonferroni, `alpha / n_tests` | family-wise error | chi-square(1) |
| `significant_fdr` | Benjamini-Hochberg over `pval` | false discovery rate | chi-square(1) |
| `significant_perm` | `--permute N` genome-wide maxima | family-wise error | nothing |
| `significant_fdr_perm` | Benjamini-Hochberg over `p_empirical` | false discovery rate | exchangeability across SNPs |

The first two convert the statistic to a p-value through a chi-square(1) null. Whether
that null fits is what **`lambda_gc`** measures — the median chi-square divided by the
chi-square(1) median, which is 1 when the null is right. IBD sharing is strongly
autocorrelated along a chromosome and pairs are not independent, so `lambda_gc` is
routinely far below 1 in practice (0.07-0.26 across regions in a 250-sample *P.
falciparum* cohort). Deflation like that means the p-values are not calibrated, and
neither correction delivers the error rate it claims. The command prints a warning
outside 0.8-1.25.

Rescaling does not fix it. Standardising the statistic by its median and MAD forces
`lambda_gc` to exactly 1 by construction, which makes the diagnostic vacuous while leaving
the tail — the part that decides significance — worse, not better.

`--permute` builds the null from the data instead. Each replicate slides every pair's IBD
segments to a random circular offset, which keeps that pair's total sharing, segment count,
segment lengths and along-genome autocorrelation intact while destroying any alignment of
segments *between* pairs. The genome-wide maximum of each replicate is one draw from the
null of "no locus is shared more than chance", and the 95th percentile of those maxima is a
real 5% family-wise threshold:

```bash
plasgenomicsutils ibd_selection_statistic --matrix ibd_matrix \
  --af freqs/allele_freqs.tsv.gz --af-group freqs/group_allele_freqs.tsv.gz \
  --meta samples.tsv --group-col region \
  --permute 200 --permute-seed 0 --output ibd_selection
```

Expect the permutation threshold to land well above the other two when `lambda_gc` is
deflated, and the significant-SNP count to fall accordingly (that cohort: 193 SNPs by
Bonferroni, 353 by BH, 78 by permutation). Cost is one full scan per replicate — roughly
4 s for 28k SNPs and 30k pairs — and each group is permuted separately, so 200 replicates
is minutes per scan. Without `--permute` the permutation columns are absent and nothing
else changes.

### Calibrated p-values, and an FDR worth quoting

The same run also writes per-SNP p-values measured against that null, which is what makes
a defensible FDR possible:

| column | meaning | resolution | assumes |
|---|---|---|---|
| `p_pointwise` | how often the null at *this SNP* reached the observed value | `1/(N+1)` | nothing beyond the shift |
| `p_empirical_binned` | how often any null *in this SNP's MAF bin* reached it | `n_bins/(N·n_snps)` | nothing beyond the shift |
| `p_empirical` | how often *any* null value anywhere reached it | `1/(N·n_snps+1)` | the null is exchangeable across MAF bins |
| `q_empirical` | Benjamini-Hochberg over whichever of the two `--empirical-pool` selects | — | the above, plus BH's usual dependence condition |

All use the Phipson-Smyth `(1 + exceedances) / (1 + draws)` form, so nothing is ever
reported as p = 0: a permutation cannot evidence a p below its own resolution.

They trade resolution against assumptions, and you have to spend one to get the other.
`p_pointwise` assumes the least but bottoms out at `1/(N+1)` — 0.005 at 200 replicates —
far too coarse to correct over ~28,000 tests; read it on the top hits. `p_empirical` pools
across every SNP to buy six orders of magnitude, which is what BH needs, but that is only
legitimate if the null has the same shape in every MAF bin.

**It often does not.** The command measures it: each bin's share of null values above one
common reference, which exchangeability says should be 0.010 everywhere. On a real
250-sample cohort those rates ran 0.000 to 0.027 across bins — some twenty standard errors
apart, so a genuine difference in null shape rather than Monte Carlo noise. Pooling then
makes `p_empirical` too *small* for SNPs in the heavy-tailed bins, which is the
anti-conservative direction. The effect is a factor of two or three, so it does not touch
the top hits, but it does move SNPs sitting near the q cutoff.

`--empirical-pool bin` drops the assumption by keeping each SNP inside its own bin
(`p_empirical_binned` drives `q_empirical` instead). That is correct rather than
approximate, and `n_bins` times coarser — with the defaults, too coarse for a genome-wide
FDR unless `N` goes up by about the same factor. Both columns are always written, so the
cheapest check is to compare them and see whether anything you care about moves.

**Pick `N` for the FDR, not just the threshold.** BH multiplies the rank-1 p-value by the
number of tests, so the smallest reachable q is about `1/N` regardless of cohort size.
With `N = 20` at `q < 0.05` the column is dead on arrival — nothing can be called. Use
`N ≥ 10/fdr_alpha` (200 at the default) for an order of magnitude of headroom. The command
computes `q_empirical_floor` and warns when it is close to your cutoff.

**What this does not fix.** `pval` and `q_value` still come from the chi-square(1) and
stay miscalibrated — do not quote them as probabilities. Use `p_empirical` / `q_empirical`
when you need a number. The *ranking* was never affected either way: every step from
`z_score` to `neg_log10_p` is monotone, so a bad reference distribution mislabels the axis
without moving any SNP relative to another.

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
