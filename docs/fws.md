# Fws (within-host diversity)

`calculate_fws` computes the per-sample Fws statistic (Manske 2012) — a monoclonal
infection scores Fws ~ 1, a polyclonal one lower — from per-sample allele depths, read
from **either a VCF/BCF or a bcftools-query AD table**. It reimplements
`moimix::getFws`.

```bash
plasgenomicsutils calculate_fws --input-vcf cohort.snps.bcf --out fws.tsv
# or from an AD table (CHROM POS REF ALT then one "ref,alt" per sample):
plasgenomicsutils calculate_fws --ad-table ad.tsv --samples samples.txt --out fws.tsv
```

## Estimators

- `--estimator regression` (default) reproduces `moimix::getFws`: `Fws = 1 − β`, the
  slope of a through-origin regression of per-sample on population heterozygosity across
  10 MAF bins. Validated against moimix on public pf7 data — see
  [Fws reconciliation](fws_moimix_reconciliation.md).
- `--estimator ratio` is a simpler summed-binned-mean estimator
  (`1 − Σ mean(Hw)/Σ mean(Hs)`).

The two are not interchangeable — don't mix a threshold tuned on one with the other's
values. moimix parity uses the defaults (`--min-depth 0 --min-alt-samples 0`).

## When Fws is not enough

Fws says how clonal a sample is, not whether one that fails the gate can still be used. An
infection whose dominant clone holds most of the parasitaemia can be re-genotyped to that clone
and treated as monoclonal; two strains of comparable size cannot. `wsaf_profile` reads that off
the allele fractions and reports, per sample, the `filter_ad_regenotype --min-freq` that would
reduce it to one clone — see [Within-host mixtures](within-host-mixtures.md).

## Notes

- `moimix::getFws` uses every biallelic record's `AD` regardless of allele string, so
  the VCF reader does too; pass a SNP-filtered callset (or `--snps-only`) for SNPs only.
- `--population-name` tags every row for later cross-cohort merging.
- `--exclude-call-regions` drops CNV windows whose within-sample heterozygosity would
  otherwise depress Fws.
