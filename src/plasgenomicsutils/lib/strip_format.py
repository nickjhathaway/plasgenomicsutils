"""Strip stale genotype-linked FORMAT fields whose stored length no longer matches the
record's genotypes.

Callers that force a genotype ploidy after calling at a different ploidy (e.g. a diploid
GT written over hexaploid calls) leave ``Number=G`` fields like ``PL`` with a value count
that disagrees with the genotypes. bcftools then refuses to subset them (``--trim-alt-alleles``
aborts with "Unexpected number of values in FORMAT/PL ... expected G=3, but found 7"), and
the values are meaningless anyway. This nulls the offending fields.

Two modes:

* ``"mismatch"`` (default) -- surgical: only records where a field's length is inconsistent
  with the per-sample genotypes are touched (the values are set missing at the correct
  length); consistent records keep their real values. Uses pysam.
* ``"always"`` -- drop the listed fields from every record and the header, via
  ``bcftools annotate -x`` (fast, streaming).

Defaults to ``PL``. Handles ``Number=G`` (from ploidy), ``Number=R``/``Number=A`` (from the
allele count), and fixed counts; fields declared ``Number=.`` are left alone.
"""

from __future__ import annotations

import os
import tempfile
from math import comb

from .bcftools import format_tags, out_flag, q, sh

# genotype-linked FORMAT fields that a re-genotype / forced-ploidy step makes stale
GENOTYPE_LINKED_FORMAT = ("PL", "GL", "GP", "PP", "PG")


def _expected_len(number, n_alleles: int, ploidy: int):
    if number == "G":
        return comb(n_alleles - 1 + ploidy, ploidy)
    if number == "R":
        return n_alleles
    if number == "A":
        return max(n_alleles - 1, 0)
    if isinstance(number, int):
        return number
    return None  # Number=. or unknown -> cannot validate


def strip_stale_format(inp: str, out: str, *, fields=("PL",),
                       mode: str = "mismatch") -> int:
    """Strip stale genotype-linked FORMAT `fields` (default ``("PL",)``).

    ``mode="mismatch"`` nulls a field only on records where its length disagrees with the
    genotypes (values set missing at the correct length; valid records untouched).
    ``mode="always"`` drops the listed fields entirely (header + records). Returns the
    number of records modified (``-1`` for ``mode="always"``, where the field is dropped
    wholesale).
    """
    if mode not in ("mismatch", "always"):
        raise ValueError("mode must be 'mismatch' or 'always'")

    if mode == "always":
        present = [f for f in fields if f in format_tags(inp)]
        fmt = out_flag(out)
        if not present:
            sh(f"bcftools view {q(inp)} -O{fmt} -o {q(out)}", tools=("bcftools",))
            return 0
        drop = ",".join(f"FORMAT/{f}" for f in present)
        sh(f"bcftools annotate -x {drop} {q(inp)} -O{fmt} -o {q(out)}", tools=("bcftools",))
        return -1

    from pysam import VariantFile  # optional dep, imported lazily

    # pysam preserves value counts correctly in VCF text, but a BCF can retain the old
    # (stale) record vector length, so write VCF and let bcftools normalise the widths into
    # the requested output format.
    tmp = tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False).name
    modified = 0
    try:
        vin = VariantFile(str(inp))
        numbers = {f: vin.header.formats[f].number for f in fields if f in vin.header.formats}
        vout = VariantFile(tmp, "wz", header=vin.header)
        try:
            for rec in vin:
                n_alleles = len(rec.alleles)
                touched = False
                for f, number in numbers.items():
                    exp = {}
                    mismatch = False
                    for name, s in rec.samples.items():
                        gt = [a for a in (s.get("GT") or ()) if a is not None]
                        exp[name] = _expected_len(number, n_alleles, len(gt) or 2)
                        val = s.get(f)
                        if val is not None and exp[name] is not None and len(val) != exp[name]:
                            mismatch = True
                    if not mismatch:
                        continue
                    for name, s in rec.samples.items():
                        if s.get(f) is not None and exp[name] is not None:
                            s[f] = (None,) * exp[name]   # drop stale values, keep valid length
                    touched = True
                if touched:
                    modified += 1
                vout.write(rec)
        finally:
            vout.close()
            vin.close()
        sh(f"bcftools view {q(tmp)} -O{out_flag(out)} -o {q(out)}", tools=("bcftools",))
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass
    return modified
