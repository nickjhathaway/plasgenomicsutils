"""Shared test helpers."""

from __future__ import annotations

import os

import pytest


@pytest.fixture(scope="session")
def bgzip():
    """Compress a plain file in place to ``<path>.gz``, returning the new path.

    Uses pysam rather than the ``bgzip`` binary. bgzip ships with htslib, which
    `brew install bcftools` pulls in but `apt-get install bcftools` does not (Debian keeps
    it in `tabix`), so shelling out passed on macOS and failed on Linux. pysam is already a
    hard dependency and carries htslib with it, so this needs nothing on PATH.
    """
    import pysam

    def _bgzip(path) -> str:
        src = str(path)
        dest = src + ".gz"
        pysam.tabix_compress(src, dest, force=True)
        os.unlink(src)
        return dest

    return _bgzip
