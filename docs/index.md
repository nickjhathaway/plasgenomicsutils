# plasgenomicsutils

Post-processing utilities for **Plasmodium** genomics data: VCF/BCF filtering and
cross-cohort harmonization, IBD post-analysis, and the Fws within-host diversity
statistic. It is the compute-side companion to
[plasgenomicsutilsR](https://nickjhathaway.github.io/plasgenomicsutilsR/), which
visualizes the outputs.

Reference-genome facts are namespaced by species (a small reference registry), so the
algorithms generalize beyond *Plasmodium falciparum* — pass `--reference`/`--species`
to select an assembly (default `pf3d7`).

## Install

```bash
pip install git+https://github.com/nickjhathaway/plasgenomicsutils
```

Or use the dedicated conda/mamba environment (see [Installation](install.md)).

## One CLI, grouped subcommands

Everything runs through a single `plasgenomicsutils` runner with a flat, globally
unique command namespace grouped into compartments (`vcf`, `ibd`, `fws`):

```bash
plasgenomicsutils --list           # catalog of all commands
plasgenomicsutils <command> -h     # help for one command
```

See the [command reference](commands.md) for the full catalog.

## The pipeline at a glance

```
BAMs → bcftools mpileup/call → per-group BCF   (or a public/GATK core-SNP VCF)
     → QC / filtering suite (hard filters, ADS, masks, AD re-genotype, coverage, MAF)
     → harmonize_bcf → bcftools merge → biallelic clean BCF
        ├── Fws:  calculate_fws
        └── IBD:  hmmibd-rs → build_ibd_matrix → analyze_ibd_matrix
                  → compute_allele_freqs → ibd_selection_statistic
```

- [VCF filtering & harmonization](filtering.md)
- [Fws](fws.md)
- [Within-host mixtures](within-host-mixtures.md)
- [IBD post-analysis](ibd.md)
- [Allele frequencies](allele-frequencies.md)
- [Coverage QC](coverage.md)

## Methods notes

- [Strand-bias sequencing artifacts](strand_bias_artifact_exclusion.md) — detecting and
  excluding SSE fake-het calls.
- [Fws reconciliation](fws_moimix_reconciliation.md) — numerical check against
  `moimix::getFws`.
