"""Rules in a file: build a Pandera schema from YAML.

A rules file names each column, its type and its checks, so the rules for a dataset can live next
to it, be reviewed like code, and be used from the command line without writing Python::

    columns:
      smiles:
        dtype: str
        checks:
          - valid_smiles
          - allowed_elements: {elements: [H, C, N, O, F, S, Cl]}
          - heavy_atoms: {min_value: 6, max_value: 50}
          - stereo_defined: {raise_warning: true}
      pic50:
        dtype: float
        checks:
          - in_range: {min_value: 0, max_value: 14}
    checks:
      - inchikey_matches: {smiles_column: smiles, inchikey_column: inchikey}

A check is written as its name alone, or as its name with its arguments. The names are
pandera-chem's checks and a set of Pandera's own (see :data:`PANDERA_CHECKS`).
"""

from __future__ import annotations

import inspect
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import pandera.pandas as pa
import yaml

from . import activity
from . import checks as chem

# pandera-chem's checks, by the name used in a rules file
CHEM_CHECKS: dict[str, Callable[..., pa.Check]] = {
    name: getattr(chem, name)
    for name in (
        "valid_smiles",
        "allowed_elements",
        "single_fragment",
        "neutral",
        "stereo_defined",
        "heavy_atoms",
        "mol_weight",
        "no_substructures",
        "unique_molecules",
        "inchikey_format",
        "inchikey_matches",
    )
}
# the activity checks (they need pint, and say how to install it when used)
CHEM_CHECKS.update(
    {
        name: getattr(activity, name)
        for name in ("concentration_units", "pchembl_matches", "activity_in_range")
    }
)
# Pandera's own checks that a rules file may use
PANDERA_CHECKS = (
    "eq",
    "ne",
    "gt",
    "ge",
    "lt",
    "le",
    "in_range",
    "isin",
    "notin",
    "str_matches",
    "str_contains",
    "str_startswith",
    "str_endswith",
    "str_length",
)
# Pandera's own options, which every pandera-chem check passes on
PANDERA_OPTIONS = {"raise_warning", "n_failure_cases", "ignore_na", "title", "description"}
DTYPES: dict[str, Any] = {"str": str, "float": float, "int": int, "bool": bool}
COLUMN_KEYS = {"dtype", "checks", "nullable", "unique", "required", "coerce", "description"}
SCHEMA_KEYS = {"columns", "checks", "strict", "coerce", "name", "description"}


class RulesError(ValueError):
    """A rules file that cannot be turned into a schema, with the reason."""


def _check(entry: object, where: str) -> pa.Check:
    """One check from its entry in a rules file.

    Args:
        entry: A check's name, or a one-item mapping of its name to its arguments.
        where: Where the entry is, for error messages (e.g. "column 'smiles'").

    Returns:
        The check.

    Raises:
        RulesError: If the entry is malformed, names an unknown check, or has wrong arguments.
    """
    if isinstance(entry, str):
        name, args = entry, {}
    elif isinstance(entry, Mapping) and len(entry) == 1:
        name, args = next(iter(entry.items()))
        args = {} if args is None else args
        if not isinstance(args, Mapping):
            raise RulesError(f"{where}: the arguments of {name!r} must be a mapping, e.g. {{a: 1}}")
    else:
        raise RulesError(f"{where}: a check is a name, or a name with its arguments: {entry!r}")
    if name in CHEM_CHECKS:
        factory = CHEM_CHECKS[name]
    elif name in PANDERA_CHECKS:
        factory = getattr(pa.Check, name)
    else:
        known = ", ".join(sorted(CHEM_CHECKS) + sorted(PANDERA_CHECKS))
        raise RulesError(f"{where}: unknown check {name!r}. Known checks: {known}")
    if name in CHEM_CHECKS:
        allowed = set(inspect.signature(factory).parameters) - {"options"} | PANDERA_OPTIONS
        unknown = set(args) - allowed
        if unknown:
            raise RulesError(
                f"{where}: {name}: unknown arguments {sorted(unknown)}; use {sorted(allowed)}"
            )
    try:
        return factory(**dict(args))
    except (TypeError, ValueError) as error:
        raise RulesError(f"{where}: {name}: {error}") from error


def _checks(entries: object, where: str) -> list[pa.Check]:
    """The list of checks under one column, or under the table.

    Args:
        entries: The entries as read from the file (a list, or nothing).
        where: Where they are, for error messages.

    Returns:
        The checks.

    Raises:
        RulesError: If ``entries`` is not a list.
    """
    if entries is None:
        return []
    if not isinstance(entries, list):
        raise RulesError(f"{where}: 'checks' must be a list")
    return [_check(entry, where) for entry in entries]


def _column(name: str, spec: object) -> pa.Column:
    """One column of the schema from its entry in a rules file.

    Args:
        name: The column's name.
        spec: Its settings: dtype, checks, nullable, unique, required, coerce, description.

    Returns:
        The column.

    Raises:
        RulesError: If a setting is unknown or the dtype is not one of :data:`DTYPES`.
    """
    where = f"column {name!r}"
    spec = {} if spec is None else spec
    if not isinstance(spec, Mapping):
        raise RulesError(f"{where}: its settings must be a mapping")
    unknown = set(spec) - COLUMN_KEYS
    if unknown:
        raise RulesError(f"{where}: unknown settings {sorted(unknown)}; use {sorted(COLUMN_KEYS)}")
    dtype = spec.get("dtype")
    if dtype is not None and dtype not in DTYPES:
        raise RulesError(f"{where}: dtype must be one of {sorted(DTYPES)}, not {dtype!r}")
    return pa.Column(
        DTYPES.get(dtype),
        _checks(spec.get("checks"), where),
        nullable=bool(spec.get("nullable", False)),
        unique=bool(spec.get("unique", False)),
        required=bool(spec.get("required", True)),
        coerce=bool(spec.get("coerce", False)),
        description=spec.get("description"),
        name=name,
    )


def schema_from_dict(rules: Mapping[str, Any]) -> pa.DataFrameSchema:
    """Build a schema from rules already read into a dictionary.

    Args:
        rules: The rules, in the form a rules file has (see the module's description).

    Returns:
        The schema.

    Raises:
        RulesError: If the rules cannot be turned into a schema; the message says where and why.
    """
    if not isinstance(rules, Mapping):
        raise RulesError("the rules must be a mapping with a 'columns' entry")
    unknown = set(rules) - SCHEMA_KEYS
    if unknown:
        raise RulesError(f"unknown top-level entries {sorted(unknown)}; use {sorted(SCHEMA_KEYS)}")
    columns = rules.get("columns")
    if not isinstance(columns, Mapping) or not columns:
        raise RulesError("the rules need a 'columns' mapping with at least one column")
    return pa.DataFrameSchema(
        {str(name): _column(str(name), spec) for name, spec in columns.items()},
        checks=_checks(rules.get("checks"), "the table"),
        strict=bool(rules.get("strict", False)),
        coerce=bool(rules.get("coerce", False)),
        name=rules.get("name"),
        description=rules.get("description"),
    )


def load_schema(path: str | Path) -> pa.DataFrameSchema:
    """Build a schema from a YAML rules file.

    Args:
        path: The rules file.

    Returns:
        The schema.

    Raises:
        RulesError: If the file is not valid YAML or its rules cannot be turned into a schema.
        FileNotFoundError: If there is no such file.
    """
    text = Path(path).read_text(encoding="utf-8")
    try:
        rules = yaml.safe_load(text)
    except yaml.YAMLError as error:
        raise RulesError(f"{path}: not valid YAML: {error}") from error
    return schema_from_dict(rules)


EXAMPLE_RULES = """\
# Rules for a table of compounds, checked by pandera-chem.
# Each check is a name, or a name with its arguments. Run:
#   pandera-chem check compounds.csv --rules rules.yaml
columns:
  compound_id:
    dtype: str
    unique: true
  smiles:
    dtype: str
    checks:
      - valid_smiles                  # RDKit can read it
      - allowed_elements              # H B C N O F P S Cl Br I (or give your own list)
      - single_fragment               # no salts or mixtures
      - unique_molecules              # no molecule twice, however it is written
      - heavy_atoms: {min_value: 6, max_value: 50}
      - stereo_defined: {raise_warning: true}   # a warning, not a failure
  pic50:
    dtype: float
    nullable: true
    checks:
      - in_range: {min_value: 0, max_value: 14}
"""
