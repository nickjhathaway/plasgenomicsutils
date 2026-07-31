"""Reference-genome registry.

This is the ONE place species/assembly-specific facts live. IBD/Fws/filtering
math is species-independent; only reference-genome facts (chromosome lengths, a
genetic-map rate, chromosome-name conventions) differ between species. Keeping
them here — keyed by a short reference id — means adding *Plasmodium vivax*
(``pvp01``) etc. later is a data change, not a code change in the algorithms.

Every command that needs such a fact takes ``--reference`` (default ``pf3d7``)
and looks it up via :func:`get_reference`.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Reference:
    """A reference genome's post-processing-relevant facts."""

    ref_id: str
    species: str
    assembly: str
    #: core (non sub-telomeric) chromosome lengths in bp, keyed by normalised chr id
    core_chrom_lengths_bp: dict[str, int]
    #: constant genetic-map rate (bp per centiMorgan)
    bp_per_cm: float

    @property
    def full_core_bp(self) -> int:
        return int(sum(self.core_chrom_lengths_bp.values()))


# --- Plasmodium falciparum 3D7 -------------------------------------------------
# Core chromosome lengths (bp) for the 14 nuclear chromosomes of the Pf3D7
# reference assembly (chromosome ids normalised to "1".."14").
PF3D7_CORE_CHROM_LENGTHS_BP: dict[str, int] = {
    "1": 640851, "2": 947102, "3": 1067971, "4": 1200490, "5": 1343557,
    "6": 1418242, "7": 1445207, "8": 1472805, "9": 1541735, "10": 1687656,
    "11": 2038340, "12": 2271494, "13": 2925236, "14": 3291936,
}

#: Pf3D7 constant genetic-map rate (bp/cM); the long-standing P. falciparum default.
PF3D7_BP_PER_CM: float = 15000.0

PF3D7 = Reference(
    ref_id="pf3d7",
    species="Plasmodium falciparum",
    assembly="Pf3D7",
    core_chrom_lengths_bp=PF3D7_CORE_CHROM_LENGTHS_BP,
    bp_per_cm=PF3D7_BP_PER_CM,
)


_REGISTRY: dict[str, Reference] = {
    PF3D7.ref_id: PF3D7,
}

#: The default reference used when a command does not specify one.
DEFAULT_REFERENCE = "pf3d7"


def available_references() -> list[str]:
    """Sorted list of registered reference ids."""
    return sorted(_REGISTRY)


def get_reference(ref_id: str = DEFAULT_REFERENCE) -> Reference:
    """Look up a :class:`Reference` by id (case-insensitive)."""
    key = ref_id.lower()
    if key not in _REGISTRY:
        raise SystemExit(
            f"ERROR: unknown --reference '{ref_id}'. "
            f"Available: {', '.join(available_references())}"
        )
    return _REGISTRY[key]


def normalise_chr(c) -> str:
    """Normalise assorted chromosome spellings to a bare number string.

    ``'Pf3D7_07_v3'``, ``'chr7'``, ``'07'``, ``7`` -> ``'7'``. Species/assembly
    prefixes seen in Plasmodium references are stripped; this is intentionally
    permissive so mixed inputs (VCF, BED, hmmibd-rs) line up on a common key.
    """
    s = str(c)
    for prefix in ("Pf3D7_", "PvP01_", "Pf_"):
        if prefix in s:
            s = s.replace(prefix, "").split("_")[0]
            break
    s = s.replace("chr", "").lstrip("0")
    return s if s else "0"
