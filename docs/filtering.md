# VCF filtering & harmonization

Each step is a parameterized command backed by `bcftools`/`bedtools`, reporting
before/after variant counts.

```bash
plasgenomicsutils hard_qc_filter --input in.bcf --output 01.bcf            # QD/MQ/SOR/RankSums, keep PASS
plasgenomicsutils singleton_filter_add_ads --input 01.bcf --output 02.bcf  # drop singletons, add FORMAT/ADS
plasgenomicsutils tandem_repeat_mask --input 02.bcf --output 03.bcf        # --bed defaults to builtin:pf3d7_tandem_repeats
plasgenomicsutils core_region_filter  --input 03.bcf --output 04.bcf       # keep core genome
plasgenomicsutils paralog_mask        --input 04.bcf --output 05.bcf       # drop paralog/multigene families (optional, see Pipeline)
plasgenomicsutils filter_ad_regenotype --input-vcf 05.bcf --output-vcf 06.bcf  # clean low AD, re-genotype
plasgenomicsutils biallelic_snp_filter --input 06.bcf --output 07.bcf      # keep biallelic SNPs
plasgenomicsutils sample_coverage_filter --input 07.bcf --output 08.bcf
plasgenomicsutils locus_missingness_filter --input 08.bcf --output 09.bcf
plasgenomicsutils maf_filter --input 09.bcf --output 10.bcf --maf-min 0.02 --maf-max 0.98
```

Region masks (`tandem_repeat_mask`, `core_region_filter`, `paralog_mask`) take `--bed`,
a plain path or a bundled asset via `builtin:<name>` (`pf3d7_core_regions`,
`pf3d7_paralog_genes`, `pf3d7_tandem_repeats`).

## Non-variant records come first

Calling a region list reports every position in it, so a callset carries a record wherever
a sample is simply reference (`ALT` is `.`). `no_alt_filter` removes those in their own
step, before any QC rule, so `variant_counts.tsv` keeps the two reasons apart — how many
positions had nothing to call, and how many real variants failed quality:

```
  input                      437
  no_alt_filter              275     <- 162 positions were non-variant
  hard_qc_filter             175     <- 100 variants failed QC
```

**This is not cosmetic: the bias statistics are computed whether or not an ALT was
called.** A non-variant record still carries `FS`, `RPBZ`, `MQBZ` and the rest, describing
the non-reference reads that were present but not called — so a hard QC rule removes
non-variant records on its own, and without this step those removals are silently mixed in
with the failing variants. On a real 400-position region list, QC took out 41 of 162
non-variant records: 31 for strand bias, 22 for read position, a few for mapping quality.

Those are positions where non-reference evidence exists but is biased enough that no ALT
was called — arguably where a reference call is least safe. `--keep-no-alts` (or
`"params": {"keep": true}`) passes them through instead, for a fill-in workflow where
"this sample is reference here" is the answer being sought; they are counted either way.
In the default chain the end result is the same, since `biallelic_snp_filter` drops
non-variant records later regardless — what changes is where they go and what the counts
tell you.

## Hard QC on a bcftools callset

`hard_qc_filter` defaults to GATK's metrics — `QD`, `MQ`, `SOR`, `MQRankSum`,
`ReadPosRankSum` — none of which a `bcftools mpileup | bcftools call` VCF contains.
`--caller bcftools` asks the same questions of the tags bcftools writes instead:

| the question | GATK | bcftools |
|---|---|---|
| is the variant only on one strand? | `SOR > 3` | `SOR > 3`, computed from `ADF`/`ADR` |
| does the variant sit at the ends of reads? | `ReadPosRankSum < -5` | `abs(RPBZ) > 5`, `abs(SCBZ) > 5` |
| are the reads carrying it poorly mapped? | `MQRankSum < -5`, `MQ < 55` | `abs(MQBZ) > 5`, `abs(MQSBZ) > 5`, `MQ < 55` |
| is the call weak for its depth? | `QD < 20` | `QUAL/INFO/DP` (off by default — see below) |
| — | — | `abs(BQBZ)`, `MQ0F` (optional extras) |

```bash
plasgenomicsutils hard_qc_filter --caller bcftools --input in.bcf --output 01.bcf
```

As a pipeline step, `caller` is written into `--emit-default-config` at its default so it
is discoverable rather than something you have to know exists:

```json
{"name": "hard_qc_filter", "params": {"caller": "bcftools"}}
```

**The `*BZ` tags are two-sided**, and getting that backwards silently keeps the records it
should drop. GATK's rank sums are signed so that one direction is the artifact, and the
usual filter is one-sided; bcftools documents its z-scores as "closer to 0 is better", and
a variant stacked at read ends turns up as either sign, so these are tested on `abs()`.

### Strand bias is measured, not tested

bcftools writes `FS`, a Fisher p-value for strand bias — but a p-value answers *am I sure
there is a skew*, and that answer turns yes for any skew at all once enough reads are
pooled. `FS` is computed over every sample at once, so a fixed cutoff means something
different in a callset of 20 samples than in one of 400: the same site, called from the
same reads, moves toward the cutoff purely because more samples carry it.

So the strand-bias test is `SOR`, GATK's strand odds ratio, which measures **how big** the
skew is and does not move when a cohort grows. bcftools writes no `SOR`, so it is computed
from the 2×2 table in `INFO/ADF` and `INFO/ADR` and compared at the same `> 3` GATK uses —
one threshold that means the same thing in either mode.

`--strand-bias-p` adds `FS` back if you want a significance test as well. Two things to
know before relying on it: at a few hundred samples it rejects sites with no meaningful
skew at all, and `INFO/FS` is a 32-bit float, so any p below about `1e-38` is stored as
exactly `0.0` — a genuinely one-strand artifact and a merely deep site become the same
number in the file.

The same caution applies in kind, though not in degree, to the `*BZ` tags: a Mann-Whitney
z-score also grows with the reads behind it, so a `--read-pos-z` that suits a small cohort
is a stricter filter on a large one. There is no effect-size form of those tags in a
bcftools callset; if a large cohort starts failing on `RPBZ` where a small one did not,
that is the reason to check before widening the threshold.

**`QD` does not carry across.** bcftools QUAL is not on GATK's scale — a clean 40x site
called at QUAL 222 has `QUAL/DP` of 5.6, so reusing GATK's `QD < 20` would throw away a
good callset. It is therefore off by default under `--caller bcftools`; pass `--qd` to set
it on a scale you have checked. Every threshold takes `none` to switch that test off.

### Calling so the tags are there

`call_variants` runs `bcftools mpileup | bcftools call` asking for exactly what the filter
reads, so the two cannot drift apart:

```bash
plasgenomicsutils call_variants --ref Pf3D7.fasta --bam-list bams.txt \
  --regions crt_region_snps.bed --threads 8 --output crt_snps.bcf
```

**`--threads` splits the region list, it does not thread bcftools.** `bcftools mpileup
--threads` only parallelises compression; the pileup itself is one core, so the way to use
a machine is a job per chunk of regions. With `--threads 8` a 400-region list is split into
8 chunks of 50, called concurrently, then indexed and `bcftools concat`-ed back into one
file. `--chunk-size` sets a fixed size instead, which evens out uneven regions at the cost
of more jobs. Without `--regions` there is nothing to split and it runs a single job.

Splitting is not quite bit-for-bit, and it is worth knowing how. The same positions come
out, with the same genotypes, depths, allele depths and bias statistics; **QUAL can move by
a few points at a handful of records** — indels, and the odd SNP beside one — because
mpileup derives indel likelihoods and BAQ from the reads around a position, and which
neighbours share a chunk changes with the split. That is `bcftools mpileup -R` itself, not
this wrapper: cutting a region file in half by hand and concatenating reproduces it exactly.
Nothing downstream here reads QUAL, but a QUAL cutoff of your own is the one thing affected.

Calling emits **all sites in the regions**, not just variants, since a region list is
usually a list of positions to fill in and a reference call there is the answer. Pass
`--variants-only` for whole-genome calling.

Chunks keep the extension of the file they came from, because **bcftools reads the
coordinate convention off it** — a `.bed` is 0-based half-open, anything else is 1-based
`CHROM POS`. Splitting a `.bed` into extensionless pieces would silently shift every region
by one base.

**Samples are named after their files.** One BAM per sample is the usual arrangement, so
`--ignore-RG` is on by default: each alignment is one sample whatever its read groups say.
bcftools then names each sample after the *path* it was given, so the last step renames
them to the file name with `--sample-suffix` removed -- `.bam` by default:

| BAM | default | `--sample-suffix .sorted.dup.pf.bam` |
|---|---|---|
| `/tank/wgs/17017-227227.sorted.dup.pf.bam` | `17017-227227.sorted.dup.pf` | `17017-227227` |

Pass the whole trailing part you want gone. A suffix that would leave two samples with the
same name is an error, raised **before** any calling starts rather than after an hour of
it. `--no-ignore-rg` reads names from the RG `SM` tags instead, and `--no-rename-samples`
keeps the path-derived names. `--dry-run` prints the commands without running them.

The equivalent by hand, if you would rather drive it yourself — most of what the filter
reads is written by default, but `FS`, `ADF`/`ADR` and `SCR` have to be asked for:

```bash
bcftools mpileup -f REF.fa \
  -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/FS,INFO/SCR \
  IN.bam -Ou | bcftools call -m -Ob -o OUT.bcf   # add -v for variants only
```

`RPBZ`, `SCBZ`, `MQBZ`, `MQSBZ`, `BQBZ`, `MQ0F`, `MQ` and `DP` come out of `mpileup`
anyway and cannot be requested explicitly. If a callset is missing what a threshold reads,
the step **stops and names the tags** rather than running: a bcftools comparison against an
absent tag is simply false, so the filter would otherwise keep every record and report
nothing wrong.

`ADF`/`ADR` are worth having beyond this step — they are what
[`strand_bias_scan`](commands.md) reads to flag SSE fake-het artifacts per site, and
`FORMAT/SP` is the per-sample version of the same question.

### Whitelisting regions from a region filter

All three take `--keep-bed`: a BED of regions to keep whatever that filter says. Use it when a
mask is right in general but wrong somewhere specific — a few positions inside a tandem repeat
that are known to be real, or a gene the paralog list catches but you trust.

```bash
plasgenomicsutils tandem_repeat_mask --input 02.bcf --output 03.bcf \
  --keep-bed keep_these.bed
```

It is region-based, so one line covers a whole gene without enumerating its variants:

```
Pf3D7_05_v3	958000	962200	pfmdr1        # every variant in the gene survives this filter
Pf3D7_07_v3	403623	403626	pfcrt-AA76    # or just one codon
```

Three things worth knowing:

- **This exempts a variant from that one rule only.** Whitelisted variants still face every
  other filter in the chain — QC, re-genotyping, coverage, MAF. It is not a "keep no matter
  what" list.
- **The BED is 0-based half-open**, like any BED, and unlike a VCF `POS`. Getting that wrong
  shifts everything a base, so each filter reports how many variants the whitelist actually
  rescued and warns when the answer is zero — a silent no-op is otherwise indistinguishable
  from success.
- **A whitelisted region rescues whole records.** Overlapping any part of a record is enough,
  so a one-base entry keeps a variant whose `REF` is longer than one base (`bedtools` sizes a
  record as `[POS-1, POS-1+len(REF))`, which reaches past a single position).

In the pipeline config it is just another param:

```json
{ "name": "tandem_repeat_mask",
  "params": { "bed": "builtin:pf3d7_tandem_repeats", "keep_bed": "keep_these.bed" } }
```

#### Building one from amino-acid positions

Residue numbers are usually how these exceptions are known ("*pfpx1* 1701 and 1705 are fine").
`aa_intervals()` in `plasgenomicsutilsR` converts them, and its default output is already
0-based half-open — the BED convention — so it can be written straight out:

```r
cds <- read_gff_cds(paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/",
                           "gff3/plasmodium_falciparum/",
                           "Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"))
codons <- aa_intervals(data.frame(transcript_id = "pfpx1", aa_position = c(1701, 1705)), cds)
write.table(codons[, c("chrom", "start", "end", "name")], "keep_these.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

Do not pass `one_based_output = TRUE` here — that is for numbers you quote in text, and would
put every interval one base off. If a codon straddles an intron its `start`/`end` span the
intron too (`spans_intron` flags it), which for a whitelist only means exempting a little more
than the three bases; use `codon_positions` if you want exactly the coding bases.

!!! note "Step order depends on your input"
    For an already-biallelic core-SNP callset the default order works. For a raw,
    multiallelic joint callset, run `biallelic_snp_filter` **before**
    `singleton_filter_add_ads` (whose `FORMAT/ADS` and the downstream MAF assume
    biallelic sites). MAF filtering is per-population — run once per cohort.

## AD cleaning and re-genotyping

`filter_ad_regenotype` judges depth by `AD`/`ADS` (not `DP`), zeros sub-threshold
allele depths per sample, recomputes `ADS`, and re-genotypes conservatively: alleles
are ranked by depth and a heterozygote is called only when the minor allele's
within-sample frequency is ≥ `--het-min-af` (default 0.2). This preserves the
mixed-infection signal Fws/COI analyses depend on. `--restrict-to-called-alleles`
narrows the caller's existing genotypes instead of re-deriving from AD.

**Ploidy.** The default keeps the conventional diploid coding used for *Plasmodium*
(`0/1` = mixed infection). `--ploidy {1,2}` sets it explicitly and is validated against
the input ploidy per record: a value **greater** than the input errors (genotypes and
likelihoods cannot be promoted), **equal** is fine, and **less** warns and trims. Use
`--ploidy 1` for true haploid calls (the single best-supported allele).

**Stale genotype-linked fields.** Re-genotyping writes a fresh `GT`, so a caller's
`Number=G` fields (`PL`/`GL`) from a different ploidy — e.g. a diploid `GT` forced over
hexaploid calls — become inconsistent and would break
`bcftools view --trim-alt-alleles`. `filter_ad_regenotype` blanks them to a consistent
length automatically as it writes. For files that did **not** go through re-genotyping,
`strip_stale_format` does the same on its own:

```bash
plasgenomicsutils strip_stale_format --input calls.bcf --output clean.bcf            # null only the inconsistent PL records (surgical, keeps valid ones)
plasgenomicsutils strip_stale_format --input calls.bcf --output clean.bcf --mode always --fields PL GL   # drop the fields entirely
```

`biallelic_snp_filter` applies the same surgical fix before trimming, so the default
chain never chokes on these fields.

`maf_filter` takes `--maf-min`/`--maf-max`; since the bounds are usually symmetric,
`--maf-max` defaults to `1 - maf_min` when unset (so `maf_min = 0.02` gives a `[0.02, 0.98]`
window). Set both for an asymmetric window.

**Per-group MAF.** With `--meta <table.tsv> --group-col country` (a per-sample metadata
table with a `sample` column and the group column), a site is kept if its minor-allele
frequency is ≥ `--maf-min` in **any** group. It is computed on the combined VCF — the
per-group frequencies pick the *sites* to keep, then that union is applied back to the
original, so **every sample's genotypes are preserved**. A carrier whose variant is below
the threshold *within its own group* keeps its real `0/1`/`1/1` call whenever the site is
kept via another group (and a `0/0` sample stays `0/0`) — there is no split-and-merge, so
nothing is blanked. A variant that is rare in *every* group is dropped.

```bash
plasgenomicsutils maf_filter --input in.bcf --output maf.bcf \
  --maf-min 0.02 --meta samples.tsv --group-col country
```

## Pipeline

Run an ordered, config-driven chain and tally counts per step:

```bash
plasgenomicsutils filter_pipeline --emit-default-config pipeline.json
plasgenomicsutils filter_pipeline --input in.bcf --config pipeline.json --outdir filtered/
```

Each step writes `filtered/NN_<name>.bcf` (indexed) plus a `variant_counts.tsv` tally
(`step`, `kind`, `count`, `path` — `kind` separates a filter's variant count from a report's
row count and from a step that was switched off).

### Turning steps on and off

A step with `"enabled": false` stays in the config and out of the run. JSON has no comments,
so that is how an optional step stays discoverable instead of being something you have to
know exists:

```json
{ "name": "paralog_mask", "enabled": false,
  "params": { "bed": "builtin:pf3d7_paralog_genes" } }
```

**`paralog_mask` ships off.** `core_region_filter` has already removed the subtelomeric
multigene families that mismap worst; most of what `paralog_mask` would take next sits in the
core and is single-copy. Some of that genuinely misbehaves, but a lot of it is fine, and
dropping all of it costs real signal — paralogy is not by itself evidence a locus is
unusable. Set `"enabled": true` when mismapping is specifically what you are controlling for,
or run the step standalone on the region you doubt.

### Why a sample was dropped

`sample_coverage_filter` writes the coverage table its decision is read from, beside the step
that produced it — `filtered/09_sample_coverage_filter_cov_info.tsv`:

| sample | n_loci | n_covered | frac_covered | n_missing_ads | mean_ads | ads_min | frac_min | margin | dropped |
|---|---|---|---|---|---|---|---|---|---|
| RCN13010 | 1901 | 843 | 0.443 | 0 | 15.64 | 5 | 0.5 | −0.057 | True |
| RCN13048 | 1901 | 912 | 0.480 | 0 | 57.26 | 5 | 0.5 | −0.020 | True |
| RCN13065 | 1901 | 964 | 0.507 | 0 | 32.15 | 5 | 0.5 | +0.007 | False |

Sorted worst-covered first, so a surprising drop is the first thing you read. `margin` is
`frac_covered - frac_min`: negative means dropped, and near zero means the sample sat on the
threshold rather than being obviously bad — `RCN13048` above missed by 0.02 while carrying a
mean ADS of 57, i.e. deeply sequenced where it *is* covered and simply absent elsewhere.
That is the distinction `mean_ads` and `n_missing_ads` exist to make, and it usually decides
whether to lower `frac_min` or drop the sample. The same rows go to the console for anything
dropped or within 0.05 of the threshold, and `--cov-table` writes them from the standalone
command.

The table is derived from the same counts as the keep/drop decision, so it always accounts for
what happened rather than being a second measurement that might disagree.

### Keeping only monoclonal infections

`fws_filter` drops every sample whose Fws falls below a threshold, leaving the ones a single
clone dominates. It ships in the default config **switched off**, because unlike everything
else in the chain it is not a quality rule — it changes which infections the callset
describes, and that is an analysis decision:

```json
{"name": "fws_filter", "enabled": false, "params": {"fws_min": 0.95}}
```

```bash
plasgenomicsutils fws_filter --input 12_maf_filter.bcf --output monoclonal.bcf --fws-min 0.92
```

Three things about where it goes and what it touches:

**It runs at the end, not as an entry gate.** Fws reads a sample's within-host diversity
against the cohort's own allele frequencies, so it wants a callset the rest of the chain has
already filtered and re-genotyped — the AD floor in `filter_ad_regenotype` is what removes
the minor-allele noise Fws would otherwise score as a second clone.

**It drops samples and no variants.** Removing samples changes every allele frequency, so a
site that cleared a MAF or missingness bar with the whole cohort may not clear it with the
one that remains. Deciding what that means for the site set is yours, so nothing is quietly
removed here; `AC`/`AN`/`AF` are refreshed on the way out, so putting `maf_filter` and
`locus_missingness_filter` after it re-applies them to the survivors:

```json
{"name": "fws_filter",  "params": {"fws_min": 0.92}},
{"name": "maf_filter",  "params": {"maf_min": 0.02}},
{"name": "locus_missingness_filter"}
```

**A sample it cannot score is dropped, not kept.** With no usable sites there is no Fws, and
an unknown sample is not a monoclonal one — keeping it would readmit exactly what the step
exists to remove. Those are counted and named separately from the polyclonal drops.

Like `sample_coverage_filter`, it writes the table its decision came from beside the step
(`filtered/13_fws_filter_fws.tsv`, or `--fws-table` from the standalone command) with one row
per sample — `sample`, `fws`, `n_sites`, `monoclonal`, `dropped` — and prints anything
dropped or within 0.05 of the threshold to the console, so a sample that missed by a hair is
visible without opening the file. A threshold that keeps nobody is an error rather than an
empty callset, since on an unfiltered input that is the likeliest reading.

Neither this nor `sample_coverage_filter` is whitelistable: `keep_bed` names regions, and
these steps judge samples. Use `calculate_fws` when the scores themselves are what you want
rather than a filtered callset — it reports the same numbers with
`--monoclonal-threshold`.

The final callset's **SNP panel BED** is written automatically next to the last step
(`filtered/NN_<last>.snps.bed`) — this is the panel the IBD tools read, so it feeds straight
into `build_ibd_matrix`:

```bash
plasgenomicsutils build_ibd_matrix --blocks hmmibd.tsv.gz \
  --snps filtered/10_maf_filter.snps.bed --snp-format bed --output ibd_matrix.npz
```

The BED's `chrom:pos` name column matches the SNP ids you'd get from the VCF, so a panel
taken as BED or VCF is interchangeable. Pass `--no-snp-bed` to skip it.

## Harmonization

Align the allele sets of separately-called cohorts before `bcftools merge`:

```bash
plasgenomicsutils harmonize_bcf --files cohortA.bcf cohortB.bcf cohortC.bcf --stub harmonized
```

It strips stale `Number=A/R/G` fields that a reshaped `AD` would otherwise break, builds
the per-site ALT union, and re-genotypes only changed samples. Indel-context records are
dropped by default (`--keep-indels` to disable).

It reports what it did per file — records read, ALT alleles zeroed and removed, records
reduced to ref-only, then records written, records that gained ALTs, and **union sites absent
from the file**.

### Sites one cohort never called

That last count is the one to read. Harmonizing settles which *alleles* the files agree on, not
which *sites* each one contains: a site only cohort A called is simply not written for cohort B.
`bcftools merge` then fills B's samples there with a missing genotype **and a missing
`FORMAT/AD`**, whole-cohort at a time.

That breaks any tool reading `AD` as an integer. hmmibd-rs stops at the first one with:

```
Error: error in processing bcf file: BcfReaderError(NumericaValueEmptyInt)
```

To keep only the sites every cohort called:

```bash
bcftools view -e 'FMT/AD="."' merged.bcf -Ob -o merged.nomiss.bcf
```

Use `FMT/AD="."`, not `FMT/AD[0]="."` — indexing a `Number=R` field per sample does not test
for a missing value and silently matches nothing. Note this is not the same as filtering on
missing genotypes: plenty of sites carry genuine per-sample missing `GT` while `AD` is present,
so `-i 'N_PASS(GT="mis")=0'` would throw away far more than necessary.

The alternative, if you would rather keep those sites, is a caller that tolerates missing depth
— `locus_missingness_filter` will trim them by missingness fraction, but any site absent from a
whole cohort fails that too.

## Strand-bias artifacts

Illumina sequence-specific errors produce fake heterozygous calls whose ALT reads sit
almost entirely on one strand. Scan a biallelic callset carrying per-strand depths and
emit a blacklist BED (feed it to `tandem_repeat_mask --bed`):

```bash
bcftools mpileup -a FORMAT/ADF,FORMAT/ADR -f ref.fa -b bams.txt -Ou \
  | bcftools norm -m- -f ref.fa -Ob -o strand.bcf
plasgenomicsutils strand_bias_scan --input-vcf strand.bcf --out-tsv verdicts.tsv --out-bed sse.bed
```

Characterize one site at read level (and optionally dump ALT reads):

```bash
plasgenomicsutils strand_read_check --bam sample.bam --pos Pf3D7_12_v3:975431 \
  --ref ref.fa --alt-base A --extract-reads
```

See [Strand-bias artifacts](strand_bias_artifact_exclusion.md) for the detection rules.
