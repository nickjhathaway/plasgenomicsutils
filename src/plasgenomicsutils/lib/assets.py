"""Resolve bundled reference BED assets shipped with the package.

Region masks can be passed a plain path, or a bundled asset by name via the
``builtin:<name>`` scheme (e.g. ``builtin:pf3d7_core_regions``), so the common
Pf3D7 masks work out of the box while any custom BED is still accepted.
"""

from __future__ import annotations

from importlib import resources

_ASSET_FILES = {
    "pf3d7_core_regions": "pf3d7_core_regions.bed",
    "pf3d7_paralog_genes": "pf3d7_paralog_genes.bed",
    "pf3d7_tandem_repeats": "pf3d7_tandem_repeats.bed",
}


def available_assets() -> list[str]:
    """The ``builtin:`` BED assets that ship with the package, as name -> path."""
    return sorted(_ASSET_FILES)


def asset_path(name: str) -> str:
    """Absolute path to a bundled asset (accepts the bare name or ``builtin:name``)."""
    key = name[len("builtin:"):] if name.startswith("builtin:") else name
    if key not in _ASSET_FILES:
        raise SystemExit(
            f"ERROR: unknown builtin asset '{key}'. Available: {', '.join(available_assets())}"
        )
    with resources.as_file(resources.files("plasgenomicsutils.assets") / _ASSET_FILES[key]) as p:
        return str(p)


def resolve_bed(value: str) -> str:
    """Map a ``--bed`` value to a filesystem path, expanding ``builtin:`` references."""
    return asset_path(value) if value.startswith("builtin:") else value
