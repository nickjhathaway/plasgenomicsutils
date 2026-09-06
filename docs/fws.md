# Fws (within-host diversity)

`calculate_fws` computes the per-sample Fws statistic (Manske 2012) — a monoclonal
infection scores Fws ~ 1, a polyclonal one lower — from per-sample allele depths, read
from **either a VCF/BCF or a bcftools-query AD table**. It reimplements
`moimix::getFws`, generalised to multiallelic sites.

```bash
plasgenomicsutils calculate_fws --input-vcf cohort.bcf --out fws.tsv
# or from an AD table (CHROM POS REF ALT then the full "ref,alt[,alt2..]" AD per sample):
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' cohort.bcf > ad.tsv
plasgenomicsutils calculate_fws --ad-table ad.tsv --samples samples.txt --out fws.tsv
```

Pass the callset **unsplit** — see [Multiallelic sites](#multiallelic-sites).

## Which sites are scored

By default only SNPs: records whose REF and every ALT are a single base. `--no-snps-only`
scores every record that carries `AD` (indels, MNP-encoded SNPs), which is what
`moimix::getFws` does. `fws_filter` has the same default and flag.

Before a site is classified, an ALT that **no sample has reads for** is dropped — the
`bcftools view --trim-alt-alleles` idea, done on `AD` rather than `GT`. A callset
joint-called across a larger cohort (pf7, say) carries alleles none of the samples in
hand support; they contribute nothing to the estimate, since zero depth is zero
frequency, but untrimmed they would make a SNP site look like an indel site to the SNP
rule and a biallelic site look multiallelic. A record with no ALT reads at all is
dropped. The counts are reported. `--no-trim` keeps every allele in the ALT column and
is needed for moimix parity, since moimix keeps the read-less records in its first MAF
bin.

## Estimators

- `--estimator regression` (default) reproduces `moimix::getFws`: `Fws = 1 − β`, the
  slope of a through-origin regression of per-sample on population heterozygosity across
  10 MAF bins. Validated against moimix on public pf7 data — see
  [Fws reconciliation](fws_moimix_reconciliation.md).
- `--estimator ratio` is a simpler summed-binned-mean estimator
  (`1 − Σ mean(Hw)/Σ mean(Hs)`).

The two are not interchangeable — don't mix a threshold tuned on one with the other's
values. moimix parity uses the defaults (`--min-depth 0 --min-alt-samples 0`).

## Multiallelic sites

Heterozygosity is computed over **every allele** at a site, `1 − Σ p²`, both within a
sample and in the population, and sites are binned on the population minor-allele
fraction (the share of reads not on the major allele). At a biallelic site this is exactly
`2p(1−p)` and `min(p, 1−p)`, i.e. moimix; at a triallelic site it is the same quantity
written for three alleles, so a sample carrying two *different* non-reference alleles is
seen as mixed.

That is why the input has to be unsplit. `bcftools norm -m-` keeps REF plus one ALT per
split record and discards the other alleles' reads, and moimix reads only the first two
`AD` columns, so under either a 50:50 mix of two ALTs looks homozygous at both records —
the error runs in the direction that lets polyclonal samples through an Fws gate:

| sample AD (ref,alt1,alt2) | true Hw | after `norm -m-` (two records) |
|---|---|---|
| 0,10,10 | 0.50 | 0.00 and 0.00 |
| 10,5,5 | 0.625 | 0.444 and 0.444 |

`--multiallelic` chooses the treatment; the count of affected records — those with reads
on more than one ALT *in this cohort*, after trimming — is reported.

- `collapse` (default): count every allele's depth, as above.
- `skip`: drop multiallelic records. Identical to `collapse` on a biallelic-only callset,
  and what moimix effectively does on a `norm -m-` split one.

Symbolic ALTs (`<NON_REF>`, `<*>`) are never counted as alleles; their `AD` column is
dropped. The spanning-deletion placeholder `*` **is** counted: reads carrying a deletion
over the site are a distinct haplotype, so a sample split between them and a base is
mixed (moimix counts them too).

## When Fws is not enough

Fws says how clonal a sample is, not whether one that fails the gate can still be used. An
infection whose dominant clone holds most of the parasitaemia can be re-genotyped to that clone
and treated as monoclonal; two strains of comparable size cannot. `wsaf_profile` reads that off
the allele fractions and reports, per sample, the `filter_ad_regenotype --min-freq` that would
reduce it to one clone — see [Within-host mixtures](within-host-mixtures.md).

## Notes

- `--population-name` tags every row for later cross-cohort merging.
- `--exclude-call-regions` drops CNV windows whose within-sample heterozygosity would
  otherwise depress Fws.
