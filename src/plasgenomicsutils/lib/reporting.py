"""How much a run says about itself.

The per-sample lines are the most useful output there is on a cohort of fifty and the
least useful on a cohort of five thousand: a coverage table with one line per sample
buries the summary that says how many were dropped. So the enumerations are a separate
level from the summaries rather than a separate flag on each step, and the level is set
once for the run.

* ``quiet`` -- nothing. The tables and files are still written; only the narration goes.
* ``verbose`` (default) -- one line per step and the summary of every report: how many
  samples were flagged, how many dropped, the counts per step.
* ``very-verbose`` -- adds the per-sample enumerations behind those summaries.
"""

from __future__ import annotations

QUIET, VERBOSE, VERY_VERBOSE = 0, 1, 2

#: Accepted ``--verbosity`` values, in order.
LEVELS = {"quiet": QUIET, "verbose": VERBOSE, "very-verbose": VERY_VERBOSE}

_level = VERBOSE


def set_verbosity(name: str) -> None:
    """Set the run's verbosity from a ``LEVELS`` name."""
    global _level
    if name not in LEVELS:
        raise SystemExit(f"ERROR: verbosity must be one of {', '.join(LEVELS)}, not {name!r}")
    _level = LEVELS[name]


def verbosity() -> int:
    return _level


def say(msg: str = "") -> None:
    """A summary line: the step, the count, how many samples were dropped."""
    if _level >= VERBOSE:
        print(msg)


def detail(msg: str = "") -> None:
    """A per-sample (or per-record) line behind a summary."""
    if _level >= VERY_VERBOSE:
        print(msg)


def listing(names, limit: int = 0) -> str:
    """Render a name list for a summary line: in full only when asked for detail.

    At ``verbose`` the count is the message and the names are the noise, so they are left
    to the table the step has already written; naming a handful is still worth it when
    there are only a handful.
    """
    names = list(names)
    if not names:
        return ""
    if _level >= VERY_VERBOSE or len(names) <= limit:
        return ": " + ", ".join(names)
    return ""
