# Allele frequencies

`compute_allele_freqs` reads a BCF/VCF once and writes what fraction of the cohort carries
each allele. It feeds `ibd_selection_statistic`, but it is also the quickest way to see
what a callset actually holds:

```bash
# every sample in the file
plasgenomicsutils compute_allele_freqs --bcf clean.bcf --output freqs/

# and per group, when there is metadata to group by
plasgenomicsutils compute_allele_freqs --bcf clean.bcf \
  --meta samples.tsv --group-col region --output freqs/
```

Without `--meta` only `allele_freqs.tsv.gz` is written; with it,
`group_allele_freqs.tsv.gz` follows.

This is not the same number as `INFO/AF`. `AF` is whatever the caller wrote at calling
time; this is computed from the genotypes as they now stand, after filtering and
re-genotyping. On a callset that has been through `filter_ad_regenotype` the two genuinely
disagree.

## Three readings of "how common is this allele"

In a polyclonal infection they differ, so all three are reported side by side.

| columns | what they count |
|---|---|
| `af`, `maf`, `ac`, `an` | alleles, over called alleles — the hard `GT` call |
| `af_weighted`, `n_samples_ad` | the mean of each sample's *within-sample* frequency |
| `prevalence`, `n_samples_alt`, `n_samples` | samples whose **genotype** carries it |
| `prevalence_ad`, `n_samples_alt_ad` | samples whose **reads** support it |

**`af`** is the classical estimate. A genotype is a hard call, so a sample carrying an
allele at 5% within-host counts exactly like one carrying it at 95%.

**`af_weighted`** keeps that information: each sample contributes its own `FORMAT/AD`
fraction, averaged over *samples* rather than over alleles. On a four-sample site with
within-host fractions 0.95, 0.05, 0.50 and 0.00, hard calls give `af = 0.625` while
`af_weighted` gives `0.375`. The denominator is each sample's own depth, so a 500x sample
does not outvote a 50x one.

**`prevalence`** is the number usually quoted for a resistance marker, and it is not `af`:
the same site above is present in 3 of 4 samples, or 0.75.

**`prevalence_ad`** asks the reads instead of the genotype, and finds the minor clones a
caller left out. On a five-sample site where one is called hom-alt but two more carry the
allele at 3% and 5% of their reads, `prevalence` is 0.20 and `prevalence_ad` is 0.60. A
sample counts as carrying the allele at `--ad-min-reads` **and** `--ad-min-freq` (2 and
0.01, matching `filter_ad_regenotype`, so "present" means the same thing in both) — both
floors, so a single stray read does not count.

```bash
# ignore the lowest-frequency clones
plasgenomicsutils compute_allele_freqs --bcf clean.bcf --output freqs/ --ad-min-freq 0.05
```

`--no-weighted` skips the AD read for speed; without usable AD the AD-derived columns are
`NaN` with one note saying so.

A frequency without its denominator cannot be read — 1.0 from two called alleles and 1.0
from seven hundred are not the same claim — so `ac`/`an` and the sample counts travel with
them. `an` also explains a `NaN`: nothing was called there.

## Multiallelic sites

By default every ALT is collapsed into one row, and `n_alts` is the only column that says
so. Three sites with 1, 2 and 3 ALTs can otherwise read identically:

```
snp_id      n_alts    af   maf  ac  an  af_weighted  prevalence
chr1:9           1 0.625 0.375   5   8        0.375        0.75
chr1:19          2 0.625 0.375   5   8        0.615        0.75
chr1:29          3 0.750 0.250   6   8        0.675        0.75
```

`--per-alt` splits them, adding `alt` and `alt_index`:

```bash
plasgenomicsutils compute_allele_freqs --bcf clean.bcf --per-alt --output freqs/
```

The split is consistent with the collapsed form: per-ALT `ac` sums to the collapsed `ac`,
per-ALT `af_weighted` sums to the collapsed one, `an` is the site's on every row, and
`n_alts` stays the site total so a row reads as "1 of 3". A biallelic site reads the same
either way.

Leave it off for tables `ibd_selection_statistic` joins on — it expects one row per SNP.

## Which frequency the selection statistic scores against

`ibd_selection_statistic --af-col` picks the column. The default is `af`, and that is the
right default: the IBD matrix comes from hmmibd-rs in `dominant-allele` mode, which reduces
each sample to a single majority allele, so the frequency of those same hard calls is what
the expected-sharing model is about. `af_weighted` describes within-host composition and
counts minor clones the IBD analysis never saw — a different question.

On a fully monoclonal cohort the two are the same quantity: with one strain per sample each
within-sample fraction is 0 or 1, so averaging them over samples is averaging hard calls.
They separate only as polyclonality rises. The flag exists so that can be checked on the
cohort at hand rather than assumed:

```bash
plasgenomicsutils ibd_selection_statistic --matrix ibd_matrix \
  --af freqs/allele_freqs.tsv.gz --af-col af_weighted --output ibd_selection_weighted
```

Naming a column an older table does not have is an error listing the columns it does have.

## Coordinates

`snp_id` is `chr:pos0` — **0-based**, matching the IBD matrix columns — and the file
carries a `#snp_coord_system=` stamp on its first line. Readers check it and refuse a table
without one, because an off-by-one panel still joins successfully and silently pairs each
SNP with the wrong frequency.
