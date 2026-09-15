"""Per-ALT carrier counts: how many **samples** carry each alternate.

"Is this variant near-private" is a question about an allele, not about a record. The two
readings agree at a biallelic site and part company as soon as there are two alternates: a
record where each ALT is private to a different sample has two singletons and no
well-supported allele, but counting samples that are merely *not* homozygous reference sees
several non-reference samples and calls the record well supported.

Written here rather than as a ``bcftools`` expression because no single expression covers
both ploidies:

* ``AC`` counts **alleles**, so a diploidized ``1/1`` carrier contributes 2.
* ``AC_Hom/2 + AC_Het`` is the right count on diploid calls -- and both tags are **zero** on
  genuinely haploid ones, which ``filter_ad_regenotype --ploidy 1`` emits.

Counting the samples directly sidesteps the question. A sample is a carrier of allele *k* if
its genotype names *k* anywhere, so a ``1/2`` het is one carrier of ALT1 and one of ALT2 --
not two carriers of either.
"""

from __future__ import annotations

#: INFO tag written by :func:`add_alt_sample_counts`: ``Number=A``, one count per ALT.
ALT_SAMPLE_TAG = "AC_SAMP"

#: Companion scalar: the best-supported **real** allele's carrier count. This is what a
#: singleton filter has to test, and ``MAX(AC_SAMP)`` is not it.
#:
#: ``*`` is a spanning deletion -- it says the sequence is not there on that haplotype. It is
#: a real observation, and it is why ``spanning_del_filter`` is default-off, but it is not a
#: variant allele. On a record like ``REF=A ALT=T,*`` where T is carried by one sample and
#: ``*`` by thirty, ``MAX(AC_SAMP)`` is 30 and the record sails through a filter whose whole
#: purpose is to drop records with no well-supported alternate. The SNP it claims to hold is
#: a singleton.
#:
#: So this tag takes the maximum over the alternates that are not ``*``. A record whose ONLY
#: alternate is ``*`` has no real allele to fall back on, and there the deletion *is* the
#: variant, so the star's own count stands -- otherwise a well-attested deletion would be
#: dropped as a singleton, which is a different rule than the one being applied.
ALT_SAMPLE_MAX_TAG = "AC_SAMP_MAX"

_SPANNING_DEL = "*"

_HEADER = {
    "ID": ALT_SAMPLE_TAG, "Number": "A", "Type": "Integer",
    "Description": "Samples carrying each ALT allele (a het counts once for each allele "
                   "it names); ploidy-independent, unlike AC",
}

_MAX_HEADER = {
    "ID": ALT_SAMPLE_MAX_TAG, "Number": "1", "Type": "Integer",
    "Description": "Largest AC_SAMP over the non-* alternates (the * count itself when * is "
                   "the only alternate); what a per-allele singleton filter tests",
}


def add_alt_sample_counts(inp: str, out: str) -> dict[str, int]:
    """Add ``INFO/AC_SAMP`` -- one carrier count per ALT -- and return a small tally.

    Records with no ALT get no tag: there is nothing to count, and a ``Number=A`` tag of
    length zero is not writable. Existing values are overwritten, so this is safe to re-run
    after anything that changes the genotypes.

    Written with pysam rather than cyvcf2 because cyvcf2's INFO setter takes only scalars and
    strings: handing it a list raises, and handing it a comma-joined string stores the tag as
    **text** despite the ``Type=Integer`` header. It reads back correctly through
    ``bcftools query``, and then ``MAX()`` fails on it with ``Unexpected type 7`` -- a
    silent-looking mismatch that only shows up at the point of use.
    """
    import pysam

    with pysam.VariantFile(inp) as vin:
        if ALT_SAMPLE_TAG not in vin.header.info:
            vin.header.info.add(_HEADER["ID"], _HEADER["Number"], _HEADER["Type"],
                                _HEADER["Description"])
        if ALT_SAMPLE_MAX_TAG not in vin.header.info:
            vin.header.info.add(_MAX_HEADER["ID"], _MAX_HEADER["Number"],
                                _MAX_HEADER["Type"], _MAX_HEADER["Description"])
        st = {"records": 0, "records_with_alt": 0, "records_all_alts_private": 0,
              "records_with_star": 0, "records_star_only_support": 0}
        with pysam.VariantFile(out, "w", header=vin.header) as vout:
            for rec in vin:
                st["records"] += 1
                n_alt = len(rec.alts or ())
                if not n_alt:
                    vout.write(rec)
                    continue
                st["records_with_alt"] += 1
                counts = [0] * n_alt
                for s in rec.samples.values():
                    named = {a for a in (s.get("GT") or ()) if a is not None and a > 0}
                    for a in named:
                        if a <= n_alt:
                            counts[a - 1] += 1
                alts = tuple(rec.alts or ())
                real = [c for a, c in zip(alts, counts) if a != _SPANNING_DEL]
                if len(real) != n_alt:
                    st["records_with_star"] += 1
                best = max(real) if real else max(counts)
                if max(counts) <= 1:
                    st["records_all_alts_private"] += 1
                if best <= 1 < max(counts):
                    # the case the scalar exists for: something is well supported here, and
                    # it is not one of the record's real alternates
                    st["records_star_only_support"] += 1
                rec.info[ALT_SAMPLE_TAG] = tuple(counts)
                rec.info[ALT_SAMPLE_MAX_TAG] = best
                vout.write(rec)
    return st
