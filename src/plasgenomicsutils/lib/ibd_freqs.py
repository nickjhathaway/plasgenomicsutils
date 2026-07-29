"""Global + per-region alternate-allele frequencies from a BCF/VCF.

Alt AF = alt-allele-count / non-missing-allele-count, per ``chr:pos``. The global
table and every region table are computed in a single pass over the file,
accumulating per-region counts as it goes. Region table is region-major with SNP
order following record order.
"""

from __future__ import annotations

import pandas as pd


def compute_allele_freqs(
    bcf_path: str,
    sample_to_region: dict[str, str] | None = None,
    zero_based: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Single pass over ``bcf_path``.

    Parameters
    ----------
    sample_to_region:
        Mapping of sample name -> region for samples present in the metadata.
        Samples absent from the mapping contribute to the global AF only (they
        are excluded from every region). ``None`` computes global AF only.
    zero_based:
        If True, emit ``chr:pos-1`` snp_ids (0-based) to match matrix labels
        from :mod:`plasgenomicsutils.lib.ibd_matrix`.

    Returns
    -------
    (global_df, region_df) where global_df has columns [snp_id, af] and region_df
    has columns [region, snp_id, af]. region_df is empty if no mapping is given.
    """
    import pysam  # lazy: only the allele-freq command needs it

    vcf = pysam.VariantFile(bcf_path)
    bcf_samples = list(vcf.header.samples)

    regions: list[str] = []
    region_of_sample: dict[str, str] = {}
    if sample_to_region:
        region_of_sample = {s: sample_to_region[s] for s in bcf_samples if s in sample_to_region}
        regions = sorted(set(region_of_sample.values()))

    global_rows: list[dict] = []
    # region -> list of {region, snp_id, af} rows, kept in record order
    region_rows: dict[str, list[dict]] = {r: [] for r in regions}

    for rec in vcf.fetch():
        pos = rec.pos - 1 if zero_based else rec.pos  # pysam pos is 1-based
        snp_id = f"{rec.chrom}:{pos}"

        g_ac = g_an = 0
        r_ac = {r: 0 for r in regions}
        r_an = {r: 0 for r in regions}

        for sname, sample in rec.samples.items():
            gt = sample["GT"]
            reg = region_of_sample.get(sname)
            for allele in gt:
                if allele is None:
                    continue
                g_an += 1
                is_alt = allele > 0
                if is_alt:
                    g_ac += 1
                if reg is not None:
                    r_an[reg] += 1
                    if is_alt:
                        r_ac[reg] += 1

        global_rows.append({"snp_id": snp_id, "af": g_ac / g_an if g_an else float("nan")})
        for r in regions:
            an = r_an[r]
            region_rows[r].append({
                "region": r, "snp_id": snp_id,
                "af": r_ac[r] / an if an else float("nan"),
            })

    vcf.close()

    global_df = pd.DataFrame(global_rows)
    if regions:
        region_df = pd.concat([pd.DataFrame(region_rows[r]) for r in regions], ignore_index=True)
    else:
        region_df = pd.DataFrame(columns=["region", "snp_id", "af"])
    return global_df, region_df
