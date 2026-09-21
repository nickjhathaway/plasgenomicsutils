"""Rewrite GATK's balanced indel pairs as the SNPs they are.

HaplotypeCaller assembles a haplotype and then aligns it to the reference with Smith-Waterman
scored at match 200, mismatch -150, gap open -260, gap extend -11. Under those weights a
multi-base substitution is often cheaper to write as an insertion and a deletion a few
bases apart than as the mismatches it is: the *pfcrt* CVIET haplotype -- ATG AAT AAA to
ATT GAA ACA over codons 74-76 -- comes out as ``403618 A>AT``, ``403622 AT>A`` and
``403625 A>C``, where three mismatches score 150 and the two gaps score 458. Nothing in
GATK then asks whether the gaps cancel, so the callset carries a substitution as two indel
records, and an SNP-only chain removes the codon 74 and 75 changes while keeping K76T.

GATK does record that the three belong together: every carrier has the same ``PID`` (and
``PGT`` phase) on all three, the physical phasing from the assembled haplotype. That is what
makes the rewrite safe. Records sharing a ``PID`` are a cluster; for each sample the phased
alleles are applied to the reference over the cluster's span, and if every haplotype comes
out the same length as the reference, the cluster is a substitution block: it is replaced by
one SNP record per mismatched base, with each sample's genotype read off its haplotypes. A
cluster where some haplotype has a net length change is a genuine indel and is left alone.

In a 249-sample sWGA cohort this is not a corner case: 920 such blocks in the core outside
tandem repeats, hiding 2,492 SNPs.

A cluster is only rewritten when at most ``max_unresolved`` samples with a call cannot
be laid out. What the derived records carry. ``AD``/``DP`` per sample come from the nearest source record
-- the cluster's records were called from one haplotype assembly and carry near-identical
depths -- with the sample's alt reads assigned to the ALT its haplotype carries. ``QUAL`` is
the minimum over the cluster and INFO is copied from the nearest source record, so the QC
metrics a later filter reads are the ones GATK computed for that haplotype. ``PL``, ``GQ``,
``PGT``, ``PID``, ``PS`` and ``SB`` are dropped: they describe genotype likelihoods of records
that no longer exist. ``INFO/SHIFTED`` names the source records. AC/AN/AF are refilled.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from collections import defaultdict

from .bcftools import index_vcf, out_flag, q, require, sh
from .reporting import say

#: FORMAT fields kept on a derived SNP record; everything else is genotype-likelihood or
#: phasing information about records that no longer exist.
KEEP_FORMAT = ("GT", "AD", "DP")
DROP_FORMAT = ("PL", "GQ", "PGT", "PID", "PS", "SB")


def _cluster_by_pid(records):
    """Group a chromosome's records by shared PID. Records with no PID are singletons."""
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    pid_first = {}
    for i, r in enumerate(records):
        parent.setdefault(i, i)
        for s in r.samples.values():
            pid = s.get("PID")
            if not pid or pid == ".":                     # pysam hands back "." as a string
                continue
            if pid in pid_first:
                union(pid_first[pid], i)
            else:
                pid_first[pid] = i
    groups = defaultdict(list)
    for i in range(len(records)):
        groups[find(i)].append(i)
    return [sorted(g) for g in groups.values() if len(g) > 1]


def _split_by_gap(idx, records, max_gap):
    """Break a PID group where consecutive records lie further apart than ``max_gap``.

    A PID names the first variant of an assembled haplotype, and assembly windows are
    hundreds of bases wide, so one PID can tie together records that have nothing to do
    with each other as a substitution block. Splitting on gaps keeps a cluster to the
    records that could share a shifted alignment.
    """
    out, cur = [], [idx[0]]
    for a, b in zip(idx, idx[1:]):
        if records[b].pos - (records[a].pos + len(records[a].ref)) > max_gap:
            out.append(cur)
            cur = [b]
        else:
            cur.append(b)
    out.append(cur)
    return [g for g in out if len(g) > 1]


def _haplotypes(sample_calls, recs, ref_seq, start):
    """A sample's two haplotype sequences over the cluster, or None if unresolvable.

    ``sample_calls`` is ``[(gt, pgt)]`` per record. Homozygous calls need no phase; a
    heterozygous call at more than one record needs PGT on every heterozygous record, or
    which alt sits with which is unknown and the sample is set missing on the derived SNPs.
    """
    haps = [[], []]
    cur = start
    n_het = sum(1 for gt, _ in sample_calls if gt and None not in gt and len(set(gt)) > 1)
    for (gt, pgt), r in zip(sample_calls, recs):
        if gt is None or None in gt:
            return None
        if r.pos < cur:                                  # overlapping records: give up
            return None
        gap = ref_seq[cur - start:r.pos - start]
        if len(set(gt)) == 1:
            alleles = (gt[0], gt[0])
        else:
            if n_het > 1:
                if not pgt or "|" not in pgt:
                    return None
                a, b = pgt.split("|")
                alleles = (int(a), int(b))
            else:
                alleles = tuple(sorted(gt))
        for h, al in enumerate(alleles):
            haps[h].append(gap + r.alleles[al])
        cur = r.pos + len(r.ref)
    tail = ref_seq[cur - start:]
    return ["".join(h) + tail for h in haps]


def _header_reference(vcf) -> str | None:
    for line in str(vcf.header).splitlines():
        if line.startswith("##reference="):
            ref = line.split("=", 1)[1].strip()
            return ref[len("file://"):] if ref.startswith("file://") else ref
    return None


def resolve_shifted_indels(inp: str, out: str, *, reference: str | None = None,
                           max_span: int = 50, max_unresolved: int = 1) -> dict:
    """Replace balanced indel clusters with the SNPs they encode. Returns counts.

    ``reference`` is the FASTA the callset was called against, needed for the bases between
    a cluster's records; when None the header's ``##reference`` line is used. ``max_span``
    bounds a cluster's reference footprint; wider clusters are left alone, since the
    reconstruction assumes one assembled haplotype per sample across the cluster, which
    GATK's physical phasing guarantees only locally.

    ``max_unresolved`` is how many samples the rewrite may cost. A sample heterozygous at
    two or more of a cluster's records with no ``PGT`` to phase them cannot be laid out, and
    writing it missing would take a genotype from a sample that has reads there; a cluster
    that would cost more than this many such samples is left exactly as it was, its records
    staying as indels. Only samples with a real call are counted -- one already no-call at
    the source records costs nothing.

    The trade is very uneven, which is why the default is 1 rather than 0 or unbounded. On a
    249-sample cohort, of 491 resolvable clusters: 195 cost nothing and give 1,495 SNPs, the
    next 60 cost one sample each and give 458 more, and the 87 worst cost over 25 samples
    each for 332 SNPs between them. At 1 the rewrite gives 1,953 SNPs for 60 genotypes out
    of ~122,000. Set 0 to never write a missing genotype at all.

    A callset without ``FORMAT/PID`` (not GATK's, or phasing stripped) is copied through
    with a note: there is nothing to cluster on.
    """
    import pysam

    require("bcftools")
    vcf = pysam.VariantFile(inp)
    counts = dict(clusters=0, resolved=0, snps_written=0, records_removed=0,
                  skipped_net_indel=0, skipped_star=0, skipped_span=0, skipped_no_change=0,
                  skipped_unresolved=0, samples_unresolved=0, samples_kept_by_skipping=0)
    if "PID" not in vcf.header.formats:
        sh(f"bcftools view {q(inp)} -O{out_flag(out)} -o {q(out)}", tools=("bcftools",))
        say("NOTE: no FORMAT/PID in the header (not a GATK callset, or phasing was stripped); "
            "nothing to resolve")
        return counts
    ref_path = reference or _header_reference(vcf)
    if not ref_path or not os.path.exists(ref_path):
        raise SystemExit(
            "resolve_shifted_indels: a reference FASTA is needed for the bases between a "
            "cluster's records; the header names none that exists"
            + (f" ({ref_path})" if ref_path else "") + ". Pass --reference.")
    return _resolve(vcf, inp, out, ref_path, max_span, max_unresolved, counts)


def _resolve(vcf, inp, out, ref_path, max_span, max_unresolved, counts):
    import pysam

    fa = pysam.FastaFile(ref_path)
    hdr = vcf.header.copy()
    hdr.add_line('##INFO=<ID=SHIFTED,Number=1,Type=String,Description="Derived from GATK '
                 'balanced-indel records at these positions (resolve_shifted_indels)">')
    tmp = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False).name
    writer = pysam.VariantFile(tmp, "w", header=hdr)
    samples = list(vcf.header.samples)

    def flush(chrom_records):
        removed = set()
        derived = []
        for group in _cluster_by_pid(chrom_records):
          for idx in _split_by_gap(group, chrom_records, max_span):
                recs = [chrom_records[i] for i in idx]
                counts["clusters"] += 1
                if not any(len(r.ref) != len(a) for r in recs for a in (r.alts or ())):
                    continue                                  # SNP-only phase group
                # A record whose only ALT is `*` is a placeholder for a deletion in this
                # cluster: it goes with the deletion, and takes no part in the haplotype.
                # A `*` beside real alternates is a genuine overlap, left alone.
                placeholders = [i for i, r in zip(idx, recs) if r.alts == ("*",)]
                recs = [r for r in recs if r.alts != ("*",)]
                if len(recs) < 2 or any(a == "*" for r in recs for a in (r.alts or ())):
                    counts["skipped_star"] += 1
                    continue
                if not any(len(r.ref) != len(a) for r in recs for a in (r.alts or ())):
                    continue
                start = recs[0].pos
                end = max(r.pos + len(r.ref) for r in recs)
                if end - start > max_span:
                    counts["skipped_span"] += 1
                    continue
                chrom = recs[0].chrom
                ref_seq = fa.fetch(chrom, start - 1, end - 1)
                haps = {}
                net_indel = False
                for s in samples:
                    calls = [(r.samples[s].get("GT"), r.samples[s].get("PGT")) for r in recs]
                    h = _haplotypes(calls, recs, ref_seq, start)
                    if h is None:
                        haps[s] = None
                        continue
                    if any(len(x) != len(ref_seq) for x in h):
                        net_indel = True
                        break
                    haps[s] = h
                if net_indel:
                    counts["skipped_net_indel"] += 1
                    continue
                # every resolved haplotype is reference-length: a substitution block
                positions = sorted({i for h in haps.values() if h for x in h
                                    for i, (a, b) in enumerate(zip(ref_seq, x)) if a != b})
                if not positions:
                    counts["skipped_no_change"] += 1
                    continue
                # A sample heterozygous at two or more of the cluster's records, with no
                # PGT to say which alternate sits on which haplotype, cannot be laid out.
                # Writing it missing would take a genotype away from a sample that has
                # reads there, so by default the whole cluster is left as it was: the
                # records stay, as indels, exactly as if this step had not run. Raising
                # `max_unresolved` trades that for the SNPs -- on a 249-sample cohort,
                # tolerating one sample recovers 458 SNPs across 60 clusters for 60
                # genotypes, while the 87 worst clusters cost over 25 samples each and
                # yield 332 SNPs between them.
                lost = [s for s, h in haps.items() if h is None
                        and any(g.get("GT") and None not in g["GT"]
                                for g in (r.samples[s] for r in recs))]
                if len(lost) > max_unresolved:
                    counts["skipped_unresolved"] += 1
                    counts["samples_kept_by_skipping"] += len(lost)
                    continue
                counts["resolved"] += 1
                counts["samples_unresolved"] += len(lost)
                removed.update(idx)                       # placeholders included
                src_positions = ",".join(str(r.pos) for r in recs)
                for i in positions:
                    pos = start + i
                    ref_base = ref_seq[i]
                    alts = sorted({x[i] for h in haps.values() if h for x in h if x[i] != ref_base})
                    nearest = min(recs, key=lambda r: min(abs(r.pos - pos), abs(r.pos + len(r.ref) - 1 - pos)))
                    new = writer.new_record(contig=chrom, start=pos - 1, stop=pos,
                                            alleles=(ref_base, *alts), qual=min(r.qual or 0 for r in recs),
                                            filter=list(nearest.filter.keys()) or None)
                    for k, v in nearest.info.items():
                        if k in hdr.info:
                            try:
                                new.info[k] = v
                            except (TypeError, ValueError):
                                pass
                    new.info["SHIFTED"] = src_positions
                    allele_index = {a: j + 1 for j, a in enumerate(alts)}
                    for s in samples:
                        h = haps.get(s)
                        call = new.samples[s]
                        src = nearest.samples[s]
                        if h is None:
                            call["GT"] = (None, None)
                            continue
                        g = tuple(0 if x[i] == ref_base else allele_index[x[i]] for x in h)
                        call["GT"] = g
                        call.phased = True
                        ad = src.get("AD")
                        if ad is not None and ad[0] is not None:
                            alt_reads = sum(a or 0 for a in ad[1:])
                            new_ad = [ad[0]] + [0] * len(alts)
                            for gi in set(g):
                                if gi:
                                    new_ad[gi] = alt_reads
                            call["AD"] = tuple(new_ad)
                        if src.get("DP") is not None:
                            call["DP"] = src["DP"]
                    derived.append(new)
                    counts["snps_written"] += 1
        counts["records_removed"] += len(removed)
        keep = [r for i, r in enumerate(chrom_records) if i not in removed]
        for r in sorted(keep + derived, key=lambda r: r.pos):
            writer.write(r)

    buf = []
    chrom = None
    for r in vcf:
        if chrom is not None and r.chrom != chrom:
            flush(buf)
            buf = []
        chrom = r.chrom
        buf.append(r)
    if buf:
        flush(buf)
    writer.close()
    fmt = out_flag(out)
    drop = ",".join(f"FORMAT/{t}" for t in DROP_FORMAT if t in vcf.header.formats)
    strip = f"| bcftools annotate -x {drop} -Ou " if drop else ""
    sh(f"bcftools sort {q(tmp)} -Ou {strip}| bcftools +fill-tags -O{fmt} -o {q(out)} -- -t AC,AN,AF",
       tools=("bcftools",))
    os.unlink(tmp)
    say(f"NOTE: {counts['resolved']:,} balanced indel cluster(s) rewritten as "
        f"{counts['snps_written']:,} SNP record(s), replacing {counts['records_removed']:,} "
        f"record(s); left alone: {counts['skipped_net_indel']:,} with a net length change, "
        f"{counts['skipped_star']:,} with a * allele, {counts['skipped_span']:,} wider than "
        f"{max_span} bp, {counts['skipped_no_change']:,} that cancel to reference"
        + (f", {counts['skipped_unresolved']:,} that would have cost more than "
           f"{max_unresolved} called sample(s) their genotype "
           f"({counts['samples_kept_by_skipping']:,} genotype(s) kept that way)"
           if counts["skipped_unresolved"] else "")
        + (f"; {counts['samples_unresolved']:,} sample(s) set missing for lack of phase"
           if counts["samples_unresolved"] else ""))
    return counts
