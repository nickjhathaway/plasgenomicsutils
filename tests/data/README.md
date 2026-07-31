# Test fixtures

## `ghana_cambodia.pf7.tiny.bcf`  (+ `.csi`, `cambodia.samples.txt`, `ghana.samples.txt`)

A small, **public** fixture for exercising the filtering pipeline on real data.

**Provenance.** Derived from the MalariaGEN Pf7 open callset, Ghana and Cambodia
2016–2017 samples only. Downsampled deterministically to be commit-safe (~2 MB):

- **Samples:** 30 per country (60 total), evenly spaced over each country's sorted
  sample list. Split lists shipped alongside for per-country runs (MAF is per-population).
- **Variants:** every 500th record genome-wide (coordinate-sorted) → ~3,100 sites
  across all 14 chromosomes. `AC/AN/AF/F_MISSING` recomputed for the 60-sample subset.
- **FORMAT:** `GT,AD,DP,GQ,PL,SB`; **INFO** includes `QD,MQ,SOR,FS,MQRankSum,ReadPosRankSum`.

Because it is a genome-wide sample of a near-raw callset it is rich in multiallelic
and `*` spanning-deletion records; after `bcftools norm -m-` + biallelic-SNP filtering
it yields ~3,600 biallelic SNPs (~570 polymorphic per country). Region coverage spans
core + non-core (subtelomere), paralog and paralog∩core, and tandem-repeat sites — so it
exercises `hard_qc_filter`, biallelic, core, paralog, tandem, coverage, missingness and
MAF. Carrying `PL` also readies it for PL-based genotype re-derivation, and the per-country
splits (trim each to its own ALTs → divergent ALT sets) drive the harmonize tests.

**Rebuild recipe:** `<homebase>/testdata_work/build_pf7_fixture.sh`
(`N_PER_COUNTRY`/`VAR_STRIDE` overridable; reads the full Ghana/Cambodia BCF; not committed).

## Unpublished data

Any unpublished callsets used for local-only testing **must never be committed** — only
the public fixtures above belong in this directory.
