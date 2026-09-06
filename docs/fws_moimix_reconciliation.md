# Fws reconciliation against `moimix::getFws`

`plasgenomicsutils calculate_fws --estimator regression` computes the Fws within-host
diversity statistic (Manske 2012), matching `moimix::getFws`. This note records the
numerical check that it is faithful, on public data, so it can be re-run.

## Result

On the committed fixture `tests/data/ghana_cambodia.pf7.tiny.bcf` (60 Ghana/Cambodia
pf7 samples), split to biallelic with `bcftools norm -m-` (11,829 records), computing
Fws over the **same site set** (every biallelic record carrying `AD`):

| | value |
|---|---|
| samples compared | 60 |
| max \|moimix − python\| (6-decimal output) | **4.9 × 10⁻⁷** |
| max \|moimix − python\| (full precision, via the Python API) | 6.7 × 10⁻¹⁶ |

The residual in the written table is the 6-decimal output rounding; at full precision the
estimator is identical. (Re-checked 2026-09-05 after the multiallelic generalisation.)

### Why the site set has to match

`moimix::getFws` uses the first two `AD` columns (ref, alt) of **every** variant in
the GDS — it does not filter by allele string, so `bcftools norm -m-` output that
still contains MNP-encoded SNPs (e.g. `REF=CG ALT=TG`) and indels is all included.
`calculate_fws` scores SNPs only by default, so parity needs `--no-snps-only`; it also
drops ALTs no sample has reads for, and a split record whose ALT has no reads at all with
it, whereas moimix keeps such records in its first MAF bin — so parity needs `--no-trim`
too. (Trimmed, the 60-sample fixture loses 8,427 read-less split records and moves by at
most 1.8 × 10⁻⁴.) Restricting one tool but not the other is what produced the initial
spurious ~0.05 differences.

### Multiallelic sites are where the two part company

`calculate_fws` counts every allele at a multiallelic site (`1 − Σ p²`; the default,
`--multiallelic collapse`) and so wants the **unsplit** callset. moimix reads only the
first two `AD` columns, and `bcftools norm -m-` discards the other alleles' reads when
it splits, so on a split callset a sample mixing two different ALTs looks homozygous to
both. Parity therefore holds on split input, where no multiallelic record survives and
`collapse` and `skip` are the same thing; on the unsplit fixture the estimates differ by
design. See [Fws](fws.md#multiallelic-sites).

At biallelic sites the generalised quantities are computed the way moimix computes them
(`1 − (p² + q²)`; minor-allele fraction as a ratio of read counts, matching
`min(coverage / sum(coverage))`), so the biallelic path is unchanged.

## Estimator choice (important)

`moimix::getFws` computes Fws as `1 − β`, where β is the slope of a **regression
through the origin** of the per-sample binned heterozygosity means on the population
binned means (`lm(sample_het ~ pop_het - 1)`), over 10 equal MAF bins. This is
`--estimator regression` (the default).

The `--estimator ratio` mode instead uses a **ratio of sums**
(`1 − Σ mean(Hw) / Σ mean(Hs)`). That is a genuinely different estimator (the
regression weights bins by population-het²), so it gives different values — do not
compare a threshold tuned on one against values from the other.

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

# python reimplementation (same site set: every split record with AD)
plasgenomicsutils calculate_fws --input-vcf fws_recon.vcf.gz \
  --estimator regression --no-snps-only --no-trim --out fws_py.tsv

# compare columns 2 of each TSV -> max|diff| ~5e-7 (6-decimal rounding)
```
