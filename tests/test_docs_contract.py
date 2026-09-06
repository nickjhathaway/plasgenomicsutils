"""The docs have to describe the CLI that exists.

A flag renamed in the parser and left behind in a code block is the kind of thing nobody
notices until someone copies the command and it fails, so these check the two against each
other rather than trusting prose.
"""

import importlib
import re
from functools import lru_cache
from pathlib import Path

import pytest

from plasgenomicsutils.cli import REGISTRY

ROOT = Path(__file__).resolve().parents[1]
DOCS = sorted((ROOT / "docs").glob("*.md")) + [ROOT / "README.md"]
COMMANDS = {c for group in REGISTRY.values() for c in group}


@lru_cache(maxsize=None)
def _parser(cmd):
    """The command's own argparse parser, built in-process.

    Each leaf module exposes `get_parser_<command>()`; going through that rather than
    shelling out to `--help` keeps this to milliseconds instead of a subprocess per
    command.
    """
    for group in REGISTRY.values():
        if cmd in group:
            mod = importlib.import_module(group[cmd].func.__module__)
            # the registry name and the module's own name are not always the same
            # (`ibd_gene_overlap` lives in `scripts/ibd/gene_overlap.py`), so find the
            # builder rather than assuming it
            builders = [v for k, v in vars(mod).items() if k.startswith("get_parser_")]
            assert len(builders) == 1, f"{cmd}: expected one get_parser_*, got {len(builders)}"
            return builders[0]()
    raise KeyError(cmd)


@lru_cache(maxsize=None)
def _help(cmd):
    return _parser(cmd).format_help()


def _documented_invocations():
    """(file, command, flag) for every flag used with a command in a ```bash block."""
    for f in DOCS:
        for block in re.findall(r"```bash\n(.*?)```", f.read_text(), re.S):
            joined = re.sub(r"\\\s*\n", " ", block)          # unwrap line continuations
            for line in joined.splitlines():
                m = re.search(r"plasgenomicsutils\s+([a-z_]+)", line)
                if not m or m.group(1) not in COMMANDS:
                    continue
                for flag in re.findall(r"\s(--[a-z0-9-]+)", line):
                    yield f.name, m.group(1), flag


def test_every_documented_flag_exists():
    flags = {}
    bad = []
    for fname, cmd, flag in _documented_invocations():
        flags.setdefault(cmd, set(re.findall(r"--[a-z0-9-]+", _help(cmd))))
        if flag not in flags[cmd]:
            bad.append(f"{fname}: `{cmd} {flag}`")
    assert bad == [], "documented flags that do not exist: " + "; ".join(bad)


def test_the_docs_actually_exercise_the_cli():
    """Guard against the check above passing because it found nothing to check."""
    seen = list(_documented_invocations())
    assert len(seen) > 40
    assert len({cmd for _, cmd, _ in seen}) > 8


def test_every_command_is_in_the_command_reference():
    listed = (ROOT / "docs" / "commands.md").read_text()
    missing = sorted(c for c in COMMANDS if c not in listed)
    assert missing == [], f"commands absent from docs/commands.md: {missing}"


@pytest.mark.parametrize("cmd", sorted(COMMANDS))
def test_every_command_has_help_and_a_description(cmd):
    out = _help(cmd)
    assert "usage:" in out, f"{cmd} has no usage line"
    # argparse prints the description between the usage block and the options block
    body = out.split("options:")[0]
    body = re.sub(r"usage:.*?\n\n", "", body, flags=re.S)
    assert len(body.strip()) > 20, f"{cmd} has no description in its --help"


# ---- every docs page is reachable from the site nav --------------------------------
# `mkdocs build --strict` fails on a page that is in docs/ but not in `nav`, since the
# config promotes omitted files to warnings. Catching it here means a new page is caught
# when it is written rather than on the next CI run.

def _nav_files(mkdocs_yml: str) -> set[str]:
    """The *.md entries of the `nav:` block.

    Read with a regex rather than a YAML parser so the guard has no dependency of its own
    and runs wherever pytest does -- the nav is ours and is a flat list of `Title: page.md`.
    """
    lines = mkdocs_yml.splitlines()
    start = next(i for i, l in enumerate(lines) if l.rstrip() == "nav:")
    out = set()
    for line in lines[start + 1:]:
        if line.strip() and not line.startswith((" ", "\t")):
            break                                   # back to a top-level key
        m = re.search(r"([\w.-]+\.md)\s*$", line)
        if m:
            out.add(m.group(1))
    return out


def test_every_docs_page_is_in_the_mkdocs_nav():
    root = Path(__file__).resolve().parents[1]
    in_nav = _nav_files((root / "mkdocs.yml").read_text())
    on_disk = {p.name for p in (root / "docs").glob("*.md")}
    missing = sorted(on_disk - in_nav)
    assert not missing, (
        f"docs/ pages absent from mkdocs.yml nav (mkdocs --strict fails on these): "
        f"{', '.join(missing)}")
    stale = sorted(in_nav - on_disk)
    assert not stale, f"mkdocs.yml nav points at pages that do not exist: {', '.join(stale)}"


def test_every_command_is_documented_under_its_own_group():
    """The catalog's groups and the reference's sections have to stay in step: a command
    filed under one group in `--list` and another in the docs is worse than ungrouped."""
    import re

    from plasgenomicsutils.cli import REGISTRY

    text = (ROOT / "docs" / "commands.md").read_text()
    sections = {}
    for m in re.finditer(r"^## `([a-z0-9_]+)`$", text, re.M):
        nxt = re.search(r"^## ", text[m.end():], re.M)
        sections[m.group(1)] = text[m.end(): m.end() + (nxt.start() if nxt else len(text))]

    assert set(REGISTRY) <= set(sections), (
        f"groups in the catalog with no section in commands.md: "
        f"{sorted(set(REGISTRY) - set(sections))}")
    for group, commands in REGISTRY.items():
        listed = set(re.findall(r"^\| `([a-z0-9_]+)` \|", sections[group], re.M))
        assert set(commands) == listed, (
            f"[{group}] differs between the catalog and commands.md: "
            f"only in catalog {sorted(set(commands) - listed)}, "
            f"only in docs {sorted(listed - set(commands))}")


def test_the_filter_pipeline_group_lists_its_steps_in_chain_order():
    """That group is ordered by the chain rather than alphabetically, on purpose -- so the
    order is worth asserting, or a later edit will quietly 'tidy' it."""
    from plasgenomicsutils.cli import REGISTRY
    from plasgenomicsutils.lib.filter_pipeline import DEFAULT_CONFIG

    listed = list(REGISTRY["vcf_filter_pipeline"])
    assert listed[0] == "filter_pipeline"          # the runner, then what it runs
    chain = [s["name"] for s in DEFAULT_CONFIG["steps"] if not s.get("report")]
    assert [c for c in listed if c in chain] == chain
