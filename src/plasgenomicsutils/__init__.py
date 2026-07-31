"""plasgenomicsutils — utilities for post processing Plasmodium genomics data.

Functional areas (CLI "compartments"):
  * vcf  — VCF/BCF filtering, AD-based re-genotyping, cross-cohort harmonization,
           and strand-bias artifact diagnostics.
  * ibd  — identity-by-descent post-analysis from hmmibd-rs blocks: binary IBD
           matrices, per-pair/per-SNP summaries, allele frequencies, and an
           EIGENSTRAT-style IBD selection statistic.
  * fws  — within-host diversity / Fws (moimix::getFws), from a VCF or AD table.

The heavy compute lives here in Python; visualization lives in the companion R
package ``plasgenomicsutilsR``.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("plasgenomicsutils")
except PackageNotFoundError:  # running from a source tree without an install
    __version__ = "0+local"

__all__ = ["__version__"]
