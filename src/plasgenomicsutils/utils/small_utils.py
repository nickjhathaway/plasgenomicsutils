"""Small shared IO / file-handling helpers used across the CLI commands."""

from __future__ import annotations

import gzip
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Iterator


class Utils:
    """Namespace of static IO/file helpers used across the CLI leaves."""

    @staticmethod
    @contextmanager
    def smart_open_read(path: str) -> Iterator[IO[str]]:
        """Open plain or gzip-compressed text for reading."""
        fh = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
        try:
            yield fh
        finally:
            fh.close()

    @staticmethod
    @contextmanager
    def smart_open_write(path: str) -> Iterator[IO[str]]:
        """Open plain or gzip-compressed text for writing; ``STDOUT`` -> stdout."""
        if path == "STDOUT" or path == "-":
            yield sys.stdout
            return
        fh = gzip.open(path, "wt") if str(path).endswith(".gz") else open(path, "w")
        try:
            yield fh
        finally:
            fh.close()

    @staticmethod
    def resolve_delim(delim: str) -> str:
        """Map the friendly ``tab``/``comma`` tokens (or a literal char) to a char."""
        if delim == "tab":
            return "\t"
        if delim == "comma":
            return ","
        return delim

    @staticmethod
    def read_table(path: str, sep=None):
        """Read a delimited metadata table into a DataFrame. When ``sep`` is ``None``,
        the delimiter is auto-detected from the header (tab if any tab is present, else
        comma), so a ``.tsv`` metadata file works without a separator flag."""
        import pandas as pd

        if sep is None:
            with Utils.smart_open_read(path) as fh:
                header = fh.readline()
            sep = "\t" if "\t" in header else ","
        return pd.read_csv(path, sep=sep)

    @staticmethod
    def output_file_check(path: str, overwrite: bool) -> None:
        """Raise unless ``path`` is writable (missing, or overwrite allowed)."""
        if path in ("STDOUT", "-"):
            return
        if Path(path).exists() and not overwrite:
            raise SystemExit(
                f"ERROR: output '{path}' already exists; pass --overwrite to replace it."
            )

    @staticmethod
    def ensure_dir(path: str) -> Path:
        """Create ``path`` (and parents) as a directory if needed; return it."""
        p = Path(path)
        p.mkdir(parents=True, exist_ok=True)
        return p

    @staticmethod
    def write_tsv_gz(df, path: str, header_comment: str | None = None) -> None:
        """Write a pandas DataFrame as a tab-delimited, gzip-compressed file.

        ``header_comment`` is written as a leading ``#`` line, for stamping metadata a
        reader should verify rather than infer (read it back with ``comment='#'``).
        """
        with gzip.open(path, "wt") as fh:
            if header_comment:
                fh.write(f"#{header_comment}\n")
            df.to_csv(fh, sep="\t", index=False)
