"""The bundled reference BED assets are shipped with the package and resolvable at
runtime from the installed location (regression test for the wheel force-include
double-add and for install-location-agnostic path resolution)."""

import os

import pytest

from plasgenomicsutils.lib.assets import asset_path, available_assets, resolve_bed

EXPECTED = {"pf3d7_core_regions", "pf3d7_paralog_genes", "pf3d7_tandem_repeats"}


def test_available_assets():
    assert set(available_assets()) == EXPECTED


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_asset_resolves_to_existing_nonempty_file(name):
    p = asset_path(name)
    assert os.path.isfile(p), p
    assert os.path.getsize(p) > 0


def test_builtin_prefix_matches_bare_name():
    assert os.path.samefile(asset_path("pf3d7_core_regions"),
                            asset_path("builtin:pf3d7_core_regions"))


def test_resolve_bed_passthrough_and_expansion():
    # a plain path is returned unchanged; a builtin: reference is expanded to a real file
    assert resolve_bed("/some/custom/mask.bed") == "/some/custom/mask.bed"
    assert os.path.isfile(resolve_bed("builtin:pf3d7_paralog_genes"))


def test_unknown_builtin_errors():
    with pytest.raises(SystemExit):
        asset_path("builtin:does_not_exist")
