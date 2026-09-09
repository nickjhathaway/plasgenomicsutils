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

## Microhaplotypes (amplicon allele tables)

Because heterozygosity is written for any number of alleles, a microhaplotype locus is
just another multiallelic site: its alleles are the haplotypes seen in the cohort and its
per-sample depths are the read counts. `--allele-table` reads the long format amplicon
pipelines write (one row per sample, locus and allele; MAD4HATTER's column names are the
defaults, `--sample-col` and friends rename them):

```bash
plasgenomicsutils calculate_fws --allele-table allele_data.tsv.gz --n-bins 0 --out fws.tsv
```

Use `--n-bins 0`. Microhaplotype loci often have no allele above 50%, so the MAF bins over
[0, 0.5] stop describing them; `--n-bins 0` regresses every locus's own `Hw` on its `Hs`
instead (`Fws = 1 − Σ Hs·Hw / Σ Hs²`). It is a different estimator from the binned one,
so do not carry a threshold between the two on SNP data.

On a Ugandan MAD4HATTER cohort (8,637 samples, 239 loci) this agreed with an Fws from the
biallelic SNPs inside the same amplicons at r = 0.996, and with WGS Fws on the 371 samples
that had both at the 0.95 gate in 94% of samples, with the same monoclonal fraction (0.962).
The disagreements were informative rather than noise: WGS is more sensitive to minor
clones at a few percent, which move a read-fraction statistic very little, and
microhaplotypes resolve mixtures of closely related strains that SNPs mostly cannot.
Amplicon read fractions carry PCR and minor-allele-filtering biases that WGS `AD` does not,
so re-check a threshold rather than assume it transfers.

`--snps-only`, `--multiallelic`, `--no-trim`, `--min-alt-samples` and
`--exclude-call-regions` do not apply to allele tables: every allele present has reads by
construction, and loci have no genomic positions or reference allele.

## Population frequencies from elsewhere

By default the population allele frequencies are the input's own pooled read fractions,
so a small or unusual batch is its own reference. `--pop-freqs` supplies them instead — to
score a handful of new samples against a reference cohort, or the same population at
another time — as a TSV with one row per locus and allele (`locus`, `allele`, `freq`;
rename with `--freq-locus-col` etc.). `--write-pop-freqs` writes that file from any input,
so the workflow is: run once on the reference cohort with `--write-pop-freqs`, then on the
new samples with `--pop-freqs`.

```bash
plasgenomicsutils calculate_fws --allele-table reference.tsv.gz --n-bins 0 \
  --write-pop-freqs ref_freqs.tsv --out ref_fws.tsv
plasgenomicsutils calculate_fws --allele-table new_batch.tsv.gz --n-bins 0 \
  --pop-freqs ref_freqs.tsv --out new_fws.tsv
```

Loci are keyed by name for `--allele-table` and by `CHROM:POS` for VCF and AD-table input,
with the allele being the haplotype or the REF/ALT string. Only the population side
changes: `Hs` and the binning variable come from the supplied frequencies, the
within-sample `Hw` never needs them. The population is exactly the alleles listed for a
locus — an allele seen here but not listed has frequency 0 there, one listed but unseen
here still counts toward `Hs`. A locus with no entry is dropped and the count reported.
Frequencies are renormalised per locus, so counts work too.

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
