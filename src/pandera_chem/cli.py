"""The command line: check a file against a rules file, without writing Python.

    pandera-chem init                                   # write a starter rules.yaml
    pandera-chem check compounds.csv --rules rules.yaml # check, write the results, report

``check`` writes the valid rows, the quarantined rows (with their reasons) and a report card (HTML
and JSON) to the output folder, prints the summary, and ends with exit code 0 if the share of
quarantined rows is within ``--max-quarantined`` (default 0: none allowed), 1 if not, and 2 if the
file or the rules cannot be used. That makes it a gate for continuous integration.
"""

from __future__ import annotations

import argparse
import contextlib
import re
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd
from pandera.errors import SchemaErrors

from ._version import __version__
from .config import EXAMPLE_RULES, RulesError, load_schema
from .quarantine import split
from .report import report_html, report_json

READERS = {".csv": ",", ".tsv": "\t", ".txt": "\t"}


class InputError(ValueError):
    """A data file that cannot be read, with the reason."""


def read_table(path: Path, smiles_column: str = "smiles") -> pd.DataFrame:
    """Read a data file into a table.

    CSV and TSV are read as text (so identifiers keep their leading zeros) and numbers are
    converted where the rules ask for them; Parquet needs ``pyarrow``; an SDF becomes one row per
    record, with its SMILES in ``smiles_column`` and its data fields as the other columns.

    Args:
        path: A .csv, .tsv, .txt (tab-separated), .parquet or .sdf file.
        smiles_column: For an SDF, the name of the column that gets each record's SMILES.

    Returns:
        The table, with a fresh index.

    Raises:
        InputError: If the file type is not supported or a needed library is missing.
        FileNotFoundError: If there is no such file.
    """
    if not path.exists():
        raise FileNotFoundError(f"no such file: {path}")
    suffix = path.suffix.lower()
    if suffix in READERS:
        return pd.read_csv(
            path, sep=READERS[suffix], dtype=str, keep_default_na=False, na_values=[""]
        )
    if suffix == ".parquet":
        try:
            return pd.read_parquet(path).reset_index(drop=True)
        except ImportError as error:
            raise InputError("reading Parquet needs pyarrow: pip install pyarrow") from error
    if suffix == ".sdf":
        return _read_sdf(path, smiles_column)
    raise InputError(
        f"cannot read {suffix or 'files without an extension'}; use .csv, .tsv, .parquet or .sdf"
    )


def _read_sdf(path: Path, smiles_column: str) -> pd.DataFrame:
    """Read an SDF: one row per record, its SMILES and its data fields.

    No record vanishes, and every record keeps its data fields (its ID, say), so a bad one can be
    traced. A record that reads but fails sanitisation keeps its unsanitised SMILES, which
    ``valid_smiles`` then quarantines; a record RDKit cannot read at all gets an empty SMILES,
    with its data fields taken from the file's text. Data fields are kept as text, as from a CSV.

    Args:
        path: The SDF.
        smiles_column: The name of the column for each record's SMILES.

    Returns:
        The table.
    """
    from rdkit import Chem, rdBase

    rows = []
    with rdBase.BlockLogs():
        supplier = Chem.SDMolSupplier(str(path), sanitize=False)
        for i in range(len(supplier)):
            mol = supplier[i]
            if mol is None:
                row = {smiles_column: ""}
                row.update(_sdf_fields(supplier.GetItemText(i)))
            else:
                row = {smiles_column: _smiles(mol)}
                row.update({name: mol.GetProp(name) for name in mol.GetPropNames()})
            rows.append(row)
    return pd.DataFrame(rows)


def _smiles(mol) -> str:
    """A record's SMILES: sanitised if it can be, otherwise as written.

    Args:
        mol: The molecule as read, not yet sanitised.

    Returns:
        The SMILES, or an empty string if none can be written.
    """
    from rdkit import Chem

    with contextlib.suppress(Exception):  # a sanitisation failure: keep the SMILES as written
        Chem.SanitizeMol(mol)
    try:
        return Chem.MolToSmiles(mol)
    except Exception:  # noqa: BLE001 - RDKit cannot write it: valid_smiles will report it
        return ""


def _sdf_fields(text: str) -> dict[str, str]:
    """The data fields of one SDF record, read from its text.

    Args:
        text: The record's text, from its first line to its ``$$$$``.

    Returns:
        Each field's name and value (multi-line values are joined with newlines).
    """
    fields = {}
    text = text.replace("\r\n", "\n")
    for match in re.finditer(r"^>.*?<([^>]+)>.*\n((?:.+\n?)*)", text, re.M):
        fields[match.group(1)] = match.group(2).strip("\n")
    return fields


def _coerce_numbers(table: pd.DataFrame, schema) -> pd.DataFrame:
    """Turn text columns into numbers where the rules give a numeric dtype.

    Text that is not a number becomes empty, and the column's ``nullable`` setting then decides
    whether that row is quarantined, so a value such as "n/a" is reported, not lost.

    Args:
        table: The table as read (text columns).
        schema: The schema from the rules.

    Returns:
        A copy with those columns converted.
    """
    out = table.copy()
    for name, column in schema.columns.items():
        if name in out.columns and str(column.dtype) in {"float64", "int64"}:
            out[name] = pd.to_numeric(out[name], errors="coerce")
            if str(column.dtype) == "int64" and out[name].notna().all():
                out[name] = out[name].astype("int64")
    return out


def check(args: argparse.Namespace) -> int:
    """Run ``pandera-chem check``.

    Args:
        args: The parsed command line.

    Returns:
        The exit code: 0 within the limit, 1 too many rows quarantined, 2 unusable input.
    """
    data, out = Path(args.data), Path(args.out)
    try:
        schema = load_schema(args.rules)
        table = _coerce_numbers(read_table(data, args.smiles_column), schema)
        result = split(table, schema)
    except (RulesError, InputError, FileNotFoundError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    except SchemaErrors as error:
        print(f"error: the table cannot be checked as a whole:\n{error}", file=sys.stderr)
        return 2
    out.mkdir(parents=True, exist_ok=True)
    stem = data.stem
    result.valid.to_csv(out / f"{stem}.valid.csv", index=False)
    result.quarantined.to_csv(out / f"{stem}.quarantine.csv", index=False)
    title = args.title or f"Data quality report: {data.name}"
    (out / f"{stem}.report.html").write_text(
        report_html(result, title, data.name), encoding="utf-8"
    )
    (out / f"{stem}.report.json").write_text(report_json(result, data.name), encoding="utf-8")

    report = result.report
    share = report.rows_quarantined / report.rows_in if report.rows_in else 0.0
    print(report)
    print(f"wrote {out.as_posix()}/{stem}.{{valid.csv,quarantine.csv,report.html,report.json}}")
    if share > args.max_quarantined:
        print(f"FAIL: {share:.1%} of rows quarantined (allowed: {args.max_quarantined:.1%})")
        return 1
    print(f"PASS: {share:.1%} of rows quarantined (allowed: {args.max_quarantined:.1%})")
    return 0


def init(args: argparse.Namespace) -> int:
    """Run ``pandera-chem init``: write a starter rules file.

    Args:
        args: The parsed command line.

    Returns:
        0, or 2 if the file exists and ``--force`` was not given.
    """
    path = Path(args.path)
    if path.exists() and not args.force:
        print(f"error: {path} already exists (use --force to replace it)", file=sys.stderr)
        return 2
    path.write_text(EXAMPLE_RULES, encoding="utf-8")
    print(f"wrote {path}; edit it, then run: pandera-chem check your-data.csv --rules {path}")
    return 0


def parser() -> argparse.ArgumentParser:
    """The command-line parser.

    Returns:
        The parser, with the ``check`` and ``init`` commands.
    """
    p = argparse.ArgumentParser(
        prog="pandera-chem",
        description="Check tables of molecules with RDKit and quarantine the bad rows.",
    )
    p.add_argument("--version", action="version", version=f"pandera-chem {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    c = sub.add_parser("check", help="check a data file against a rules file")
    c.add_argument("data", help="a .csv, .tsv, .parquet or .sdf file")
    c.add_argument("--rules", required=True, help="the YAML rules file")
    c.add_argument("--out", default="pandera-chem-out", help="output folder (default: %(default)s)")
    c.add_argument(
        "--max-quarantined",
        type=float,
        default=0.0,
        metavar="SHARE",
        help="largest share of rows allowed in quarantine before failing, e.g. 0.05 (default: 0)",
    )
    c.add_argument("--smiles-column", default="smiles", help="for an SDF: the SMILES column's name")
    c.add_argument("--title", help="the report's title")
    c.set_defaults(run=check)

    i = sub.add_parser("init", help="write a starter rules file")
    i.add_argument("path", nargs="?", default="rules.yaml", help="where (default: %(default)s)")
    i.add_argument("--force", action="store_true", help="replace the file if it exists")
    i.set_defaults(run=init)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command line.

    Args:
        argv: The arguments (default: those the program was started with).

    Returns:
        The exit code.
    """
    args = parser().parse_args(argv)
    return args.run(args)
