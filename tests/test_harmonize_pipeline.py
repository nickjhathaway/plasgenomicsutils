"""Integration test: harmonize two divergent cohorts derived from the shipped pf7 fixture.

Splits the fixture by country and trims each to its own ALT alleles, producing two
cohorts with divergent ALT sets (each "missing" the other's alts). Harmonize should
reconcile them to a common ALT union and drop the stale PL field so bcftools merge
succeeds. Skips cleanly where bcftools/bedtools or the fixture are unavailable.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HAVE_BCFTOOLS = shutil.which("bcftools") is not None
HAVE_BEDTOOLS = shutil.which("bedtools") is not None
FIXTURE = Path(__file__).parent / "data" / "ghana_cambodia.pf7.tiny.bcf"

pytestmark = pytest.mark.skipif(
    not (HAVE_BCFTOOLS and HAVE_BEDTOOLS and FIXTURE.exists()),
    reason="needs bcftools+bedtools and the shipped pf7 fixture",
)


def _run(cmd):
    return subprocess.run(cmd, shell=True, executable="/bin/bash",
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def _positions_with_alt(bcf):
    out = _run(f"bcftools query -f '%CHROM:%POS\\t%REF>%ALT\\n' {bcf}").stdout
    return {line.split("\t")[0]: line.split("\t")[1] for line in out.splitlines()}


def test_harmonize_divergent_cohorts_then_merge(tmp_path):
    d = tmp_path
    data = FIXTURE.parent
    # 1. split by country + trim to each cohort's own ALTs -> divergent ALT sets
    for c in ("ghana", "cambodia"):
        r = _run(f"bcftools view -S {data/f'{c}.samples.txt'} {FIXTURE} -Ou "
                 f"| bcftools view --trim-alt-alleles -Ob -o {d/f'{c}.bcf'} "
                 f"&& bcftools index -f {d/f'{c}.bcf'}")
        assert r.returncode == 0, r.stderr

    before = _positions_with_alt(d / "ghana.bcf")
    after_cam = _positions_with_alt(d / "cambodia.bcf")
    diverging = [p for p in before if p in after_cam and before[p] != after_cam[p]]
    assert len(diverging) > 50, "expected many positions with divergent ALTs to harmonize"

    # 2. harmonize the two cohorts (BCF out -> exercises the PL strip + convert path)
    r = _run(f"{sys.executable} -m plasgenomicsutils.cli harmonize_bcf "
             f"--files {d/'ghana.bcf'} {d/'cambodia.bcf'} --stub {d/'h'} --output-format b")
    assert r.returncode == 0, r.stderr
    for c in ("ghana", "cambodia"):
        assert (d / f"h_{c}.bcf").exists()
        _run(f"bcftools index -f {d/f'h_{c}.bcf'}")
        # PL must be gone (harmonize can't recompute it after reshaping alleles)
        hdr = _run(f"bcftools view -h {d/f'h_{c}.bcf'}").stdout
        assert "ID=PL," not in hdr

    # 3. harmonized cohorts share an identical ALT set at shared positions
    hg, hc = _positions_with_alt(d / "h_ghana.bcf"), _positions_with_alt(d / "h_cambodia.bcf")
    assert [p for p in hg if p in hc and hg[p] != hc[p]] == []

    # 4. the payoff: bcftools merge succeeds (used to fail on the stale PL length)
    m = _run(f"bcftools merge {d/'h_ghana.bcf'} {d/'h_cambodia.bcf'} --merge snps "
             f"-Ob -o {d/'merged.bcf'}")
    assert m.returncode == 0, m.stderr
    n = int(_run(f"bcftools view -H {d/'merged.bcf'} | wc -l").stdout.strip())
    assert n > 0
    n_samples = len(_run(f"bcftools query -l {d/'merged.bcf'}").stdout.split())
    assert n_samples == 60
