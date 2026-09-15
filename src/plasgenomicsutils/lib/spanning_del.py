"""Recode the calls that name a ``*`` allele as missing. **Opt-in, not a default.**

``*`` in ALT says a deletion called somewhere else covers this position in some samples. It
names no base and does not say *which* deletion, so as an allele it is a poor one: two
samples carrying it cannot be said to share an origin, which is the test that makes one
alternate worth separating from another.

But it is not nothing either, and that is why this is off by default. ``*`` is a **confident
observation that the sequence is absent**, not a failure to call. A site with 20 deleted
samples and 5 carrying a variant, run through this, comes out reading as though 20 samples
could not be genotyped there -- which is not what the data says. In *P. falciparum* the
distinction matters more than usual: across the dimorphic regions a deletion is frequently
the other haplotype rather than a dropout, so which samples carry it is a result.

**Turn it on when the question is about the variants of the non-deleted strains** and the
deletion itself is not the subject -- looking at 3D7-type sequence where the alternative
haplotype is simply absent, say. Then ``*`` is genuinely in the way: left alone it is scored
as a third allele, which inflates every heterozygosity-based statistic with a non-allele and
makes an ordinary SNP under someone else's deletion look multiallelic, so a ``-M2`` filter
deletes it. On the shipped Pf7 fixture, 312 post-trim records are multiallelic **only**
because of a ``*``, and recoding cuts the genuine multiallelic burden from 544 records to 232.

Whichever way that goes, the **counting** question is separate and always applies:
:func:`~plasgenomicsutils.lib.bcftools.classify_record` no longer lets a ``*`` decide a
record's class, and :func:`~plasgenomicsutils.lib.bcftools.count_spanning_del` reports
carriage orthogonally. That is reporting only and changes no data.

The cost of the recode is bounded and is reported per run: at the 369 fixture records where
``*`` sits beside real alleles, a mean 17.2% of called samples carry it, 77% of records lose
under a quarter, and none loses more than 90%. The count has to reach the missingness filters
-- a site that quietly loses half its samples is exactly what resurfaces later as an
unexplained outlier.

A call is nulled slot by slot, not wholesale. ``*/T`` means one haplotype is deleted and the
other carries ``T``; that is a real observation and it is kept, as ``./T``. Only 85.5% of the
fixture's ``*`` calls name nothing else, and those do become fully missing.

This module only rewrites genotypes. Dropping the now-uncarried ``*`` from the ALT column, and
re-laying the ``Number=R``/``A``/``G`` fields that go with it, is
``bcftools view --trim-alt-alleles``'s job -- see
:func:`~plasgenomicsutils.lib.vcf_filters.spanning_del_filter`, which composes the two.
"""

from __future__ import annotations

#: Keys of the dict every function here returns, in report order.
SPANNING_DEL_COUNTS = (
    "records", "records_with_spanning_del", "records_star_only",
    "calls_recoded", "slots_recoded", "calls_fully_missing",
)


def _blank() -> dict[str, int]:
    return dict.fromkeys(SPANNING_DEL_COUNTS, 0)


def spanning_del_to_missing(inp: str, out: str) -> dict[str, int]:
    """Null every genotype slot that names a ``*`` allele. Returns the tally.

    The ALT column is left exactly as it was, so this is separable from the trim and can be
    checked on its own. Records with no ``*`` are written through untouched.

    Counts returned: ``records`` seen; ``records_with_spanning_del``;
    ``records_star_only`` (nothing but ``*`` in ALT, so not a variant site at all);
    ``calls_recoded`` (sample-calls with at least one slot nulled); ``slots_recoded``
    (allele slots nulled, which is larger when a call names ``*`` twice); and
    ``calls_fully_missing`` (calls that had nothing left afterwards).
    """
    from cyvcf2 import VCF, Writer

    vcf = VCF(inp)
    writer = Writer(out, vcf)
    st = _blank()
    try:
        for v in vcf:
            st["records"] += 1
            alts = list(v.ALT)
            star = [i + 1 for i, a in enumerate(alts) if a == "*"]
            if not star:
                writer.write_record(v)
                continue
            st["records_with_spanning_del"] += 1
            if all(a == "*" for a in alts):
                st["records_star_only"] += 1

            starred = set(star)
            gts = v.genotypes          # [[a1, a2, ..., phased], ...]
            changed = False
            for g in gts:
                hit = [k for k in range(len(g) - 1) if g[k] in starred]
                if not hit:
                    continue
                for k in hit:
                    g[k] = -1
                st["calls_recoded"] += 1
                st["slots_recoded"] += len(hit)
                if all(a < 0 for a in g[:-1]):
                    st["calls_fully_missing"] += 1
                changed = True
            if changed:
                v.genotypes = gts
            writer.write_record(v)
    finally:
        writer.close()
        vcf.close()
    return st


def spanning_del_note(st: dict[str, int]) -> str:
    """One line saying what the recode cost, for the step report.

    Said out loud rather than left in the return value: the samples a site loses here do not
    show up as dropped records, so nothing else in the run would mention them.
    """
    if not st.get("records_with_spanning_del"):
        return "     no `*` alleles: nothing to recode"
    parts = [
        f"{st['records_with_spanning_del']:,} record(s) carried a `*` allele",
        f"{st['records_star_only']:,} of them had no other ALT",
        f"{st['calls_recoded']:,} call(s) recoded to missing"
        f" ({st['calls_fully_missing']:,} wholly)",
    ]
    return "     " + "; ".join(parts)
