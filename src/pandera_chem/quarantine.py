"""Quarantine: keep the rows that pass, set the others aside with their reasons.

Pandera either accepts a table or raises an error listing what failed. In a pipeline you usually
want neither: you want to carry on with the good rows and keep the bad ones, each with the
reason it was set aside, so someone can fix them later. :func:`split` does that, and its report
gives counts that must add up.
"""

from __future__ import annotations

import warnings
from collections import Counter
from dataclasses import dataclass, field

import pandas as pd
import pandera.pandas as pa
from pandera.errors import SchemaErrors, SchemaWarning

REASONS = "reasons"  # the column added to the quarantined rows


@dataclass(frozen=True)
class Report:
    """What happened to a table: how many rows came in, passed and were set aside, and why.

    Attributes:
        rows_in: Rows in the table given to :func:`split`.
        rows_valid: Rows that passed every check.
        rows_quarantined: Rows set aside.
        by_reason: For each reason ("column: check"), how many rows failed it. A row that fails
            two checks is counted under both, so these can add up to more than
            ``rows_quarantined``.
        warnings: The messages of checks that only warn (``raise_warning=True``).
    """

    rows_in: int
    rows_valid: int
    rows_quarantined: int
    by_reason: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    @property
    def reconciles(self) -> bool:
        """Whether every row is accounted for: rows in = valid + quarantined.

        Returns:
            True if the counts add up.
        """
        return self.rows_in == self.rows_valid + self.rows_quarantined

    def __str__(self) -> str:
        """A short summary to print or log.

        Returns:
            The summary, one fact per line.
        """
        lines = [
            f"{self.rows_in} rows in: {self.rows_valid} valid, {self.rows_quarantined} quarantined",
        ]
        lines += [f"  {n:>6}  {reason}" for reason, n in self.by_reason.items()]
        lines += [f"  warning: {w}" for w in self.warnings]
        return "\n".join(lines)


@dataclass(frozen=True)
class Split:
    """The result of :func:`split`.

    Attributes:
        valid: The rows that passed, as validated (and coerced) by the schema.
        quarantined: The rows that failed, unchanged, with a ``reasons`` column added.
        report: The counts.
    """

    valid: pd.DataFrame
    quarantined: pd.DataFrame
    report: Report


def _reason(context: object, column: object, check: object) -> str:
    """One failure's reason: "column: check", or just the check for a table-wide one.

    Args:
        context: Where Pandera ran the check: "DataFrameSchema" for a table-wide check.
        column: The column that failed (Pandera repeats a table-wide failure for each column).
        check: The check's name.

    Returns:
        The reason.
    """
    if context == "DataFrameSchema" or pd.isna(column):
        return str(check)
    return f"{column}: {check}"


def _warnings(caught: list[warnings.WarningMessage]) -> list[str]:
    """The distinct messages of the checks that only warn, in the order they came.

    Args:
        caught: Warnings recorded while validating.

    Returns:
        Each Pandera warning's message, once.
    """
    messages = [str(w.message) for w in caught if issubclass(w.category, SchemaWarning)]
    return list(dict.fromkeys(messages))


def _all_valid(df: pd.DataFrame, valid: pd.DataFrame, caught: list) -> Split:
    """The split of a table in which every row passed.

    Args:
        df: The table given to :func:`split`.
        valid: The table as validated by the schema.
        caught: Warnings recorded while validating.

    Returns:
        Every row valid, none quarantined.
    """
    empty = df.iloc[0:0].assign(**{REASONS: pd.Series(dtype=str)})
    return Split(valid, empty, Report(len(df), len(df), 0, {}, _warnings(caught)))


def split(df: pd.DataFrame, schema: pa.DataFrameSchema) -> Split:
    """Validate a table and divide it into the rows that pass and the rows that fail.

    Args:
        df: The table. Its index must be unique, so each failure can be traced to one row.
        schema: The Pandera schema to check it against.

    Returns:
        The valid rows, the quarantined rows with their reasons, and a report.

    Raises:
        ValueError: If the index of ``df`` has repeated labels.
        pandera.errors.SchemaErrors: If something fails that no row can be blamed for, such as
            a missing column; the whole table is then unusable, and the error says why.
    """
    if not df.index.is_unique:
        raise ValueError("the table's index has repeated labels; reset it first")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SchemaWarning)
        try:
            return _all_valid(df, schema.validate(df, lazy=True), caught)
        except SchemaErrors as error:
            failures = error.failure_cases
            if failures["index"].isna().any():
                raise  # a failure of the whole table, not of particular rows
        reasons = (
            failures.assign(
                reason=[
                    _reason(s, c, k)
                    for s, c, k in zip(
                        failures["schema_context"],
                        failures["column"],
                        failures["check"],
                        strict=True,
                    )
                ]
            )
            .drop_duplicates(["index", "reason"])
            .groupby("index", sort=False)["reason"]
            .agg("; ".join)
        )
        bad = df.index.isin(reasons.index)
        quarantined = df[bad].assign(**{REASONS: df.index[bad].map(reasons)})
        valid = schema.validate(df[~bad], lazy=True) if (~bad).any() else df[~bad]
    warned = _warnings(caught)
    counts = Counter(r for joined in reasons for r in joined.split("; "))
    report = Report(len(df), int((~bad).sum()), int(bad.sum()), dict(counts.most_common()), warned)
    return Split(valid, quarantined, report)
