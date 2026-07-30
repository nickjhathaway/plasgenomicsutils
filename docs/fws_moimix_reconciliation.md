# Fws reconciliation against `moimix::getFws`

`plasgenomicsutils calculate_fws --estimator regression` reimplements the Fws
within-host diversity statistic (Manske 2012) so the R/`moimix` VCF→GDS pipeline is
no longer needed. This note records the numerical check that the reimplementation is
faithful, on public data, so it can be re-run.

## Result

On the committed fixture `tests/data/ghana_cambodia.pf7.tiny.bcf` (60 Ghana/Cambodia
pf7 samples), split to biallelic with `bcftools norm -m-` (11,829 records), computing
Fws over the **same site set** (every biallelic record carrying `AD`):

| | value |
|---|---|
| samples compared | 60 |
| max \|moimix − python\| | **4.8 × 10⁻⁵** |
| mean \|diff\| | 2.5 × 10⁻⁵ |

The residual is entirely the 4-decimal rounding of the earlier output format
(`0.9289` vs moimix `0.928911`); at full precision the estimator is identical. The
output now prints 6 decimals.

### Why the site set has to match

`moimix::getFws` uses the first two `AD` columns (ref, alt) of **every** variant in
the GDS — it does not filter by allele string, so `bcftools norm -m-` output that
still contains MNP-encoded SNPs (e.g. `REF=CG ALT=TG`) and indels is all included.
`read_ad_vcf` therefore also reads every biallelic record with `AD` (not just
single-base REF/ALT); pass a SNP-filtered VCF if you want SNPs only (or
`--snps-only`). Restricting one tool but not the other is what produced the initial
spurious ~0.05 differences.

## Estimator choice (important)

`moimix::getFws` computes Fws as `1 − β`, where β is the slope of a **regression
through the origin** of the per-sample binned heterozygosity means on the population
binned means (`lm(sample_het ~ pop_het - 1)`), over 10 equal MAF bins. This is
`--estimator regression` (the default).

The `wgs_cnv_workflow` monoclonality gate (`src/pf_cnv/fws.py`) instead uses a
**ratio of sums** (`1 − Σ mean(Hw) / Σ mean(Hs)`). That is a genuinely different
estimator (the regression weights bins by population-het²), so it gives different
values — do not compare a threshold tuned on one against values from the other. It is
preserved as `--estimator ratio`.

## Reproduce

Requires R with `moimix`, `SeqArray`, `SeqVarTools` (here: homebrew R 4.6.1), plus
`bcftools` and the `plasgenomicsutils` env.

```bash
FIX=tests/data/ghana_cambodia.pf7.tiny.bcf
bcftools norm -m- "$FIX" -Oz -o fws_recon.vcf.gz && bcftools index -t fws_recon.vcf.gz

# moimix reference
Rscript -e '
  suppressPackageStartupMessages({library(moimix); library(SeqArray)})
  seqVCF2GDS("fws_recon.vcf.gz","fws_recon.gds",verbose=FALSE)
  h<-seqOpen("fws_recon.gds"); f<-getFws(h); seqClose(h)
  write.table(data.frame(sample=names(f),fws=as.numeric(f)),
              "fws_moimix.tsv",sep="\t",quote=FALSE,row.names=FALSE)'

# python reimplementation (same site set)
plasgenomicsutils calculate_fws --input-vcf fws_recon.vcf.gz \
  --estimator regression --out fws_py.tsv

# compare columns 2 of each TSV -> max|diff| ~5e-5 (rounding)
```
