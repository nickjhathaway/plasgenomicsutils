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

## Development

```bash
git clone https://github.com/nickjhathaway/plasgenomicsutils
cd plasgenomicsutils
pip install -e ".[dev]"
pytest
```

## License

GPL-3.
