# Installation

## pip

```bash
pip install git+https://github.com/nickjhathaway/plasgenomicsutils
```

This exposes the `plasgenomicsutils` CLI. External tools `bcftools` and `bedtools`
must also be on `PATH` (several commands shell out to them).

## conda / mamba environment

A dedicated environment file pins every dependency (Python libraries plus
`bcftools`/`bedtools`/`htslib`) so runs and tests are reproducible:

```bash
mamba env create -f environment.yml
mamba activate plasgenomicsutils
pip install -e .
```

## Dependencies

- Python ≥ 3.10 with `numpy`, `pandas`, `scipy`, `pysam`, `cyvcf2`
- `bcftools`, `bedtools`, `htslib` (external CLIs)

`cyvcf2` powers the fast bulk-FORMAT paths (e.g. `filter_ad_regenotype`,
`strand_bias_scan`); `pysam` handles allele-reshaping steps such as `harmonize_bcf`.
