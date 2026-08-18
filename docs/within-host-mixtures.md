# Within-host mixtures: which polyclonal samples can still be used

[Fws](fws.md) says how clonal a sample is. It does not say whether a sample that *fails* the
gate can still be used, and that is a separate question with a concrete answer: an infection
whose dominant clone carries most of the parasitaemia can be re-genotyped to that clone and
treated as monoclonal, while two strains of comparable size cannot.

`wsaf_profile` answers the operational version — **for a sample failing the Fws gate, is there
a dominant clone, and what filter would reduce it to that clone?**

```bash
plasgenomicsutils wsaf_profile --input-vcf cohort.bcf --out wsaf.tsv \
  --sites-out wsaf_sites.tsv.gz --fws results/fws.tsv
```

Pass `--fws` and the report is focused on the samples that fail the gate, which is where the
question matters. Without it, every sample is reported the same way.

## The threshold is a statement about composition, not a tuning knob

A dominant clone holding **70%** means the minor alleles sit at **30% or below**, which is
exactly `filter_ad_regenotype --min-freq 0.30`. The two numbers are the same choice, so
`--min-dominant 0.70` is a claim about what parasite composition you are willing to call one
clone — not a parameter to tune until the counts look right.

It also makes the prediction checkable. Filtering at `f` zeroes every minor allele below `f`,
so the sites left heterozygous afterwards are exactly those whose minor fraction is at or above
`f`. Their rate over the covered sites is `residual_het_rate`, and `min_freq_needed` inverts
it: the smallest threshold that gets that residue under `--max-residual-het`.

## The columns that matter

| column | what it says |
|---|---|
| `min_freq_needed` | smallest `filter_ad_regenotype --min-freq` that reduces this sample to one clone, or `NaN` if none below 0.5 does. **The argument to pass, per sample** |
| `dominant_frac` | `1 - min_freq_needed`: what share the dominant clone holds once the rest is removed |
| `residual_het_rate` | fraction of covered sites that would stay heterozygous at the implied threshold |
| `class` | `monoclonal` / `dominant_clone` / `mixed` / `undetermined` |
| `minor_mode`, `n_bands` | the fullest band of the minor fraction, and how many bands there are |
| `minor_median`, `minor_q95` | quantiles — descriptive only, see below |
| `n_sites`, `n_het`, `het_rate` | how much of the genome is covered, and heterozygous at all |
| `wsmaf_mean` | mean within-sample frequency of the population-level minor allele, the quantity the COI literature reports |

### The classes

- **`monoclonal`** — already under the residue bar with no filtering. Not the same as no
  heterozygous calls: a deeply sequenced clonal sample accumulates a scattering of them from
  error alone, which is why a *rate* is asked for rather than a count.
- **`dominant_clone`** — a filter at or below `1 - min_dominant` reduces it to its dominant
  clone. **These are the ones worth rescuing.**
- **`mixed`** — it would take a harder filter than that, or none below 0.5 works at all. A
  minor-allele filter has to stay under 0.5 to keep the dominant call, so forcing one would
  delete a real strain and leave a chimera of the rest.
- **`undetermined`** — too few covered sites to judge. Usually thin coverage.

### Rates, not shares

Every decision is a rate over the sites a sample **covers**, never a share of its heterozygous
sites. A share carries no weight: a sample with 250 het sites of which 13% are near 0.5 has 33
such sites in a 20,000-site panel — a few repetitive or mismapped loci, not a second genome.
Two genomes in one host differ across the genome.

The quantiles are reported because they describe the distribution, but they *are* shares of het
sites, so on a sample with few of them a handful of bad loci move them a long way: a clonal
sample with 142 het calls can show `minor_q95` near 0.45 off seven junk sites. Don't threshold
on them.

## What `min_freq_needed` depends on

It rises with **both** the minor strain's proportion and how much of the genome is
heterozygous, because a filter has to clear the binomial spread around the strain proportion
rather than just its centre. Two strains differing at 40% of sites leave more residue at a given
threshold than two differing at 15% — so a strain at 25% can still need a filter above 0.30 when
the pair is heavily divergent.

## Then filter and re-gate

```bash
plasgenomicsutils filter_ad_regenotype --input-vcf cohort.bcf \
  --output-vcf cohort.dominant.bcf --min-freq 0.30
plasgenomicsutils calculate_fws --input-vcf cohort.dominant.bcf --out fws_filtered.tsv
```

Raising `--min-freq` zeroes within-sample minor alleles below that fraction, so an infection
with a smaller second clone re-genotypes as the clone that dominates it — the same reduction
`hmmibd-rs --bcf-read-mode dominant-allele` makes when it reads a BCF. The profile predicts the
outcome; the re-run confirms it, so gate on the re-run.

Two things to know:

- **Raising `--min-freq` only removes minor alleles**, so what it can do to Fws is push a sample
  toward its dominant strain. Check the before-and-after on your own cohort rather than assuming
  a threshold transfers.
- **A sample can move non-monotonically** across thresholds under the default regression
  estimator, because changing which sites are heterozygous changes which MAF bin they land in,
  which changes the regression the estimator fits. `--estimator ratio` is monotonic in the
  threshold.

## The callset has to be filtered first

These rates are only meaningful on a quality-filtered callset. On a raw one — subtelomeric and
hypervariable regions included, or multiallelic records carrying `*` spanning-deletion alleles
inside their `AD` — mismapping alone produces heterozygous calls across the genome at a rate
that swamps any real signal, and nearly every sample looks polyclonal. Biallelic SNPs in the
core genome is the intended input.

## Plotting it

`--sites-out` writes the per-`(sample, site)` fractions behind the summary, which
`plasgenomicsutilsR::plot_wsaf()` draws as one small histogram per sample, ordered and coloured
by class:

```r
plot_wsaf(readr::read_tsv("wsaf_sites.tsv.gz"), profile = readr::read_tsv("wsaf.tsv"))
```

The shape is easier to judge than any summary column: mass squeezed against zero is one
dominant strain, a band near 0.5 is a co-dominant one, and more than one band is more than two
strains. Judge by how much *area* sits past the line, not by whether anything does.

## References

- Zhu, S. J. *et al.* (2019) The origins and relatedness structure of mixed infections vary
  with local prevalence of *P. falciparum* malaria. *eLife* 8, e40845.
- Paschalidis, A. *et al.* (2023) coiaf: directly estimating complexity of infection with
  allele frequencies. *PLOS Computational Biology* 19(6), e1010247.
