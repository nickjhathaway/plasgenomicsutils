# Excluding strand-bias sequencing artifacts from VAF-based CNV estimation

**Scope:** These are sWGA / NextSeq (Illumina, 2-channel) whole-genome BAMs of monoclonal
*Plasmodium falciparum*. A downstream analysis scans a region for heterozygous-looking calls
and uses their variant allele fraction (VAF) to estimate copy number (CN). This document
describes a class of **sequencing artifact that produces fake het calls at fractional VAF**
and must be excluded, plus concrete rules to detect and drop them.

---

## 1. TL;DR — what to exclude and why

A subset of positions produce a "variant" whose alt reads are **almost entirely on one strand**
(e.g. ~45% VAF on the reverse strand, ~0% on the forward strand, pooling to ~0.25–0.33 overall).
These are **not real alleles** — they are Illumina **sequence-specific errors (SSE)** driven by
local sequence context (GC islands, G-homopolymers, GGC motifs; Nakamura et al. 2011, *NAR* 39:e90).

They are dangerous for CN-from-VAF specifically because a strand-biased artifact lands at
**exactly the fractional VAF (0.25, 0.33) that a real 4-copy or 3-copy locus would show**, so it
mimics the CNV signal it is meant to measure. Any position where the alt is strand-restricted
and/or low-quality should be excluded from VAF/CN estimation.

**Core rule:** a real het/CNV allele is **strand-balanced** (similar VAF on forward and reverse
reads). An artifact is **strand-restricted**. Compute per-strand VAF and drop the strand-restricted ones.

---

## 2. The artifact class (background)

- **Mechanism:** certain motifs — inverted repeats, GGC sequences, and G-homopolymer runs —
  cause lagging-strand **dephasing** during sequencing-by-synthesis. The cluster loses phase
  synchrony, producing miscalled bases with **low base quality**. Because the trigger is a
  physical structure encountered in one synthesis direction, the error is **strand-specific**.
- **In the AT-rich *P. falciparum* genome (~81% AT, ~19% GC), GC islands and G-runs are rare and
  stand out**, so SSE concentrates at a small number of predictable, reproducible positions.
- **Signature at the read level:**
  - alt reads almost exclusively on **one strand** (usually reverse),
  - **low base quality** on the alt base (often median BQ < 15),
  - all three non-reference bases may appear on that strand (A dominant, C/T minor) — i.e. the
    reads are *garbage-called*, not carrying a specific chemical lesion,
  - alt base spread across sequencing cycles (position-anchored, not a single cycle),
  - reproducible across samples that cover the locus in the same orientation.
- **Why it can be library/panel-specific:** the SSE substrate is genomic (present in every
  sample), but it only becomes *visible* where a library has deep, oriented coverage reading
  *into* the motif. A primer set / prep that covers the locus that way exposes it; others that
  don't cover it there never show it. So "only seen in one library" does **not** mean it is real —
  it means only that library reads into the motif.

---

## 3. Why this corrupts CN-from-VAF

CN estimation from het VAF assumes VAF ≈ (copies carrying alt) / (total copies):

| apparent VAF | naive CN interpretation |
|---|---|
| 0.50 | 2 copies, 1 alt |
| 0.33 | 3 copies, 1 alt |
| 0.25 | 4 copies, 1 alt |

A strand-restricted artifact with ~45% alt on one strand and ~0% on the other pools to
**~0.22–0.30 overall**, which is indistinguishable from a genuine 3–4 copy signal **if you only
look at pooled VAF**. The strand decomposition is what separates them. Therefore: **compute VAF
per strand, and never trust a pooled VAF whose alt support is strand-restricted.**

---

## 4. Detection signature — what to compute per candidate site (per sample)

You need **strand-resolved allele depths**. From bcftools, request them at mpileup time:

```bash
bcftools mpileup -B -d 100000 -q 0 -Q 0 \
  -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,FORMAT/DP \
  -f REF.fa -r CHROM:POS-POS --bam-list bams.txt -Ou \
| bcftools norm -m- -f REF.fa -Ou \
| bcftools view -e 'ALT="."' \
| bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%ADF\t%ADR\t%SP\n]'
```

After `norm -m-`, `ADF` and `ADR` are two-element arrays `ref,alt` per strand.

For each candidate het (per sample), compute:

- `ref_fwd, alt_fwd` = ADF; `ref_rev, alt_rev` = ADR
- `vaf_fwd = alt_fwd / (ref_fwd + alt_fwd)`
- `vaf_rev = alt_rev / (ref_rev + alt_rev)`
- `vaf     = (alt_fwd + alt_rev) / total_depth`
- **strand-bias (Fisher exact)** on the 2×2 table `[[ref_fwd, alt_fwd],[ref_rev, alt_rev]]`,
  as a Phred score `SB = -10*log10(p)`. (`FORMAT/SP` from bcftools is an equivalent Phred SB.)
- **alt-read base quality** (median BQ of alt-supporting reads) if available — a strong
  corroborating flag; requires a pileup pass (see §6).

Discriminators:

| quantity | real het / CNV allele | strand-bias artifact |
|---|---|---|
| `vaf_fwd` vs `vaf_rev` | similar (on-diagonal) | very different (one ~0) |
| minor-strand alt fraction | substantial | ~0 |
| Fisher SB (Phred) | low (< ~20) | very high (often > 60) |
| alt-read median BQ | normal (≥ 25) | low (often < 15) |

---

## 5. Exclusion rules (concrete, tunable)

Let `strand_minor = min(vaf_fwd, vaf_rev)`, `strand_major = max(vaf_fwd, vaf_rev)`.

**Primary exclusion (drop the site for this sample if ANY of):**

1. **Strand-restricted alt:** `strand_minor / strand_major < 0.15`
   *(i.e. > ~87% of the alt signal is on one strand)*, **and** the minor strand has adequate depth
   to have detected it (e.g. `ref_minorstrand + alt_minorstrand >= 20`). The depth condition
   prevents flagging genuinely low-coverage strands.
2. **Extreme Fisher strand bias:** `SB_phred > 40` (or `FORMAT/SP > 40`).
3. **Low-quality alt support:** alt-read median `BQ < 20` (when read-level BQ is available).

**Recommended combined rule (specific to the artifact, conservative about false drops):**

> Exclude if **(strand_minor/strand_major < 0.15 AND minor-strand depth ≥ 20)**
> **OR** **(SB_phred > 60)**
> **OR** **(strand_minor/strand_major < 0.30 AND alt median BQ < 15).**

Rationale: for CN-from-VAF you average VAF across many het sites in the target region, so dropping
a handful of suspicious sites costs little, while a single artifact at 0.25/0.33 can produce a
false CN call. Bias toward excluding.

**Also recommended:**

- **Maintain a BED blacklist** of confirmed artifact positions (SSE is position-reproducible), and
  hard-exclude those across all samples regardless of the per-sample numbers. Seed it with the
  worked example in §7.
- **Aggregate, don't trust singletons:** estimate CN from the distribution of VAF across all
  *surviving* (strand-balanced, high-BQ) het sites in the gene, not from any one site.

---

## 6. Implementation notes

### 6a. Per-strand VAF + Fisher SB from ADF/ADR (portable, VCF-only)

```python
from math import log10
from scipy.stats import fisher_exact

def strand_bias_verdict(ref_fwd, alt_fwd, ref_rev, alt_rev,
                        min_minor_depth=20, ratio_hard=0.15, ratio_soft=0.30,
                        sb_hard=60, alt_bq_median=None):
    dpf, dpr = ref_fwd + alt_fwd, ref_rev + alt_rev
    vaf_f = alt_fwd / dpf if dpf else 0.0
    vaf_r = alt_rev / dpr if dpr else 0.0
    minor, major = sorted((vaf_f, vaf_r))
    ratio = (minor / major) if major > 0 else 1.0
    p = fisher_exact([[ref_fwd, alt_fwd], [ref_rev, alt_rev]])[1]
    sb_phred = -10 * log10(max(p, 1e-300))
    minor_depth = dpf if vaf_f <= vaf_r else dpr

    drop, reasons = False, []
    if ratio < ratio_hard and minor_depth >= min_minor_depth:
        drop = True; reasons.append(f"strand-restricted (ratio={ratio:.2f})")
    if sb_phred > sb_hard:
        drop = True; reasons.append(f"fisher SB Phred={sb_phred:.0f}")
    if alt_bq_median is not None and ratio < ratio_soft and alt_bq_median < 15:
        drop = True; reasons.append(f"low alt BQ={alt_bq_median} + strand skew")
    return {"drop": drop, "vaf_fwd": vaf_f, "vaf_rev": vaf_r,
            "sb_phred": sb_phred, "reasons": reasons}
```

### 6b. Optional read-level base quality (pysam)

If you want the alt-read median BQ (strongest single corroborator), do a pileup pass and, for
reads whose base at POS equals the alt, collect `read.query_qualities[query_position]` and take
the median. Low median (< ~15) alongside strand skew is essentially diagnostic of SSE.

### 6c. Proactive reference-context pre-flagging (optional but cheap)

Because SSE is sequence-driven, you can flag risky positions *before* looking at reads by scanning
the reference across the target region for Nakamura triggers, which are rare in Pf and therefore
specific:

- **G/C homopolymer runs ≥ 3** (e.g. `GGGG`, `CCCC`),
- **GGC / GCC motifs**,
- **local GC content spikes** (e.g. a ≥ 20 bp window whose GC% is far above the ~19% genome-wide
  background — a "GC island").

Positions within a few bp of these, showing strand-restricted low-BQ alts, are almost certainly
SSE. This lets you build the BED blacklist deterministically rather than discovering artifacts one
at a time.

---

## 7. Worked example — the archetype to calibrate against

**Position `Pf3D7_12_v3:975431` (ref `G`), sWGA/NextSeq Ugandan cohort.** Apparent `G>A` het at
pooled VAF ~0.25–0.30. It is an artifact. Representative sample:

```
reads over site: 3716  (fwd 1465, rev 2251)
fwd ALT reads: 2      -> vaf_fwd = 0.1%
rev ALT reads: 1021   -> vaf_rev = 45.4%
pooled alt VAF ~ 27.5%              <- would falsely read as ~3-4 copies
alt-read base quality: median 9, 87% < Q20
alt base spread across cycles 4-131 (not a single cycle)
soft-clip 0%, SA-tag 0% (not a chimera)
```

- `strand_minor/strand_major = 0.001/0.454 ≈ 0.002` → far below 0.15, minor-strand depth ~1463
  (ample) → **drop**.
- Consistent across samples (also seen at ~32% and ~42% rev VAF in two others), and G>C / G>T
  appear reverse-only at lower rates too — confirming garbage basecalling, not a real allele.

**Reference context (why this position):** a GC island in an AT-rich background —
`GGGG` homopolymer at 975426–975429 (2 bp 5′ of POS) and a GC-rich patch `GAGCACG` at
975437–975443; core ±12 bp is ~44% GC vs ~19% genome-wide. Reverse reads synthesize through the
3′ GC patch and dephase onto 975431. Add `Pf3D7_12_v3:975431` to the blacklist.

---

## 8. Caveats

- **Thresholds are tunable.** The values above bias toward excluding artifacts; if you find real
  variants being dropped, relax `ratio_hard` / `sb_hard`, but keep the base-quality corroborator.
- **Low coverage ≠ strand bias.** Do not flag a site merely because one strand is shallow;
  require the minor strand to have enough depth to have *detected* the alt (the `min_minor_depth`
  guard). Genuine low-coverage strands are a coverage problem, not an SSE artifact.
- **Don't over-filter homopolymer-adjacent real variants.** Some true variants sit near G-runs;
  the strand + BQ evidence (not context alone) is what justifies exclusion. Context pre-flagging is
  a prioritization aid, not a standalone reason to drop.
- **CN from VAF should be an aggregate.** Use the surviving strand-balanced, high-BQ het sites
  across the gene; a robust summary (e.g. median VAF of clean sites) is far safer than any single
  locus.
- **Reference:** Nakamura K, et al. "Sequence-specific error profile of Illumina sequencers."
  *Nucleic Acids Research* 2011; 39(13):e90.
