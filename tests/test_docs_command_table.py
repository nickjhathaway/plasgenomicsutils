"""The command catalogue in docs/commands.md must match the CLI registry.

The table is hand-written, so it silently falls behind when a command is added -- this
pins it to the registry instead.
"""

import re
from pathlib import Path

from plasgenomicsutils.cli import REGISTRY


DOCS = Path(__file__).resolve().parents[1] / "docs" / "commands.md"


def _documented():
    text = DOCS.read_text()
    return set(re.findall(r"^\| `([a-z0-9_]+)` \|", text, flags=re.M))


def test_every_command_is_documented():
    registered = {name for cmds in REGISTRY.values() for name in cmds}
    missing = registered - _documented()
    assert not missing, f"missing from docs/commands.md: {sorted(missing)}"


def test_no_stale_commands_documented():
    registered = {name for cmds in REGISTRY.values() for name in cmds}
    stale = _documented() - registered
    assert not stale, f"documented but not registered: {sorted(stale)}"
