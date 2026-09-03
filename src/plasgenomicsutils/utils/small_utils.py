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
    def resolve_column(columns, want: str, *, source: str = "table",
                       required: bool = True):
        """The column in ``columns`` that means ``want``, matching case-insensitively.

        Metadata comes from whoever assembled it, and ``Sample`` / ``sample`` /
        ``SAMPLE`` are the same column to everyone except a string comparison. An exact
        match always wins; otherwise a single case-insensitive match is used. Two columns
        differing only in case is an error rather than a coin toss, since which one holds
        the ids is not ours to guess.
        """
        cols = list(columns)
        if want in cols:
            return want
        hits = [c for c in cols if str(c).lower() == want.lower()]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise SystemExit(
                f"ERROR: {source} has {len(hits)} columns differing only in case for "
                f"'{want}': {', '.join(map(str, hits))}. Rename all but one."
            )
        if not required:
            return None
        raise SystemExit(
            f"ERROR: column '{want}' not in {source} "
            f"(has: {', '.join(map(str, cols))})"
        )

    @staticmethod
    def normalise_columns(df, wants=("sample",), *, source: str = "metadata",
                          quiet: bool = False):
        """Rename whichever case-variant of each name in ``wants`` the table uses.

        Applied when a metadata table is read, so everything downstream can name the
        column one way. The rename is announced -- a table read differently from how it
        is written is worth one line of output.
        """
        renames = {}
        for want in wants:
            found = Utils.resolve_column(df.columns, want, source=source, required=False)
            if found is not None and found != want:
                renames[found] = want
        if renames:
            if not quiet:
                for old, new in renames.items():
                    print(f"  note: {source} column '{old}' read as '{new}'")
            df = df.rename(columns=renames)
        return df

    @staticmethod
    def normalise_sample_ids(df, col: str = "sample"):
        """Make ``df[col]`` a column of strings, if it is there.

        **A sample id is a label, not a number.** An all-numeric cohort -- micronix ids,
        barcodes, anything without a letter in it -- reads back from a metadata table as
        int64, while every id it has to be joined against (a VCF sample name, an IBD pair
        label, a matrix row) is a string. Those never match, and nothing errors: a merge
        just finds nothing, so groups come out empty and the run looks like it worked.

        Integral floats are formatted without their trailing ``.0``, since one blank cell
        is enough to make pandas read the whole column as float and turn ``4064928010``
        into ``4064928010.0``.
        """
        import pandas as pd

        if col not in df.columns or df[col].dtype == object:
            return df
        df = df.copy()
        vals = df[col]
        if pd.api.types.is_float_dtype(vals):
            df[col] = vals.map(lambda v: "" if pd.isna(v)
                               else (str(int(v)) if float(v).is_integer() else str(v)))
        else:
            df[col] = vals.astype(str)
        return df

    @staticmethod
    def read_meta(path: str, sep=None, wants=("sample",), quiet: bool = False):
        """Read a metadata table, normalise the case of its key columns, and make the
        sample ids strings -- see :func:`normalise_sample_ids` for why that matters."""
        df = Utils.normalise_columns(Utils.read_table(path, sep), wants,
                                     source=f"metadata ({path})", quiet=quiet)
        return Utils.normalise_sample_ids(df)

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
