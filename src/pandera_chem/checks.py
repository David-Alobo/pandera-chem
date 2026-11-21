"""Chemistry checks for Pandera schemas.

Each function returns a ``pandera.Check`` to put on a column of SMILES (or of InChIKeys), next
to Pandera's own checks::

    import pandera.pandas as pa
    import pandera_chem as pc

    schema = pa.DataFrameSchema({
        "smiles": pa.Column(str, [pc.valid_smiles(), pc.allowed_elements(), pc.neutral()]),
    })

One rule keeps the reasons clear: only :func:`valid_smiles` fails a SMILES that cannot be read.
Every other check passes it, so each bad record is reported once, for its real problem.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from typing import Any

import pandas as pd
import pandera.pandas as pa
from rdkit import Chem
from rdkit.Chem import Descriptors

from ._mol import inchikey, parse

# the elements of most drug-like and screening libraries; pass your own to allowed_elements()
ORGANIC = frozenset({"H", "B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"})
INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


def valid_smiles(**options: Any) -> pa.Check:
    """The SMILES can be read and sanitised by RDKit.

    Catches syntax errors, unclosed rings, impossible valences and aromatic rings that cannot
    be kekulised. :func:`pandera_chem.problems` says why a given SMILES fails.

    Args:
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """
    return pa.Check(
        lambda smiles: parse(smiles) is not None,
        element_wise=True,
        name="valid_smiles",
        description="the SMILES cannot be read or sanitised by RDKit",
        **options,
    )


def allowed_elements(elements: Iterable[str] = ORGANIC, **options: Any) -> pa.Check:
    """Every atom is one of the allowed elements.

    Args:
        elements: The element symbols allowed, e.g. ``{"C", "N", "O", "H"}``. Defaults to
            :data:`ORGANIC`: H, B, C, N, O, F, P, S, Cl, Br and I.
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """
    allowed = frozenset(elements)

    def ok(smiles: str) -> bool:
        """Whether every atom of one SMILES is allowed (unreadable SMILES pass).

        Args:
            smiles: The SMILES.

        Returns:
            True if allowed or unreadable.
        """
        mol = parse(smiles)
        return mol is None or all(a.GetSymbol() in allowed for a in mol.GetAtoms())

    return pa.Check(
        ok,
        element_wise=True,
        name="allowed_elements",
        description=f"an element outside {sorted(allowed)}",
        statistics={"elements": sorted(allowed)},
        **options,
    )


def single_fragment(**options: Any) -> pa.Check:
    """The record is one connected molecule: no salts, solvents or mixtures.

    Args:
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """

    def ok(smiles: str) -> bool:
        """Whether one SMILES is a single piece (unreadable SMILES pass).

        Args:
            smiles: The SMILES.

        Returns:
            True if one piece or unreadable.
        """
        mol = parse(smiles)
        return mol is None or len(Chem.GetMolFrags(mol)) == 1

    return pa.Check(
        ok,
        element_wise=True,
        name="single_fragment",
        description="more than one molecule in the record (a salt, solvent or mixture)",
        **options,
    )


def neutral(**options: Any) -> pa.Check:
    """The net formal charge is zero (charges that cancel, as in a nitro group, are fine).

    Args:
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """

    def ok(smiles: str) -> bool:
        """Whether one SMILES has no net charge (unreadable SMILES pass).

        Args:
            smiles: The SMILES.

        Returns:
            True if neutral or unreadable.
        """
        mol = parse(smiles)
        return mol is None or Chem.GetFormalCharge(mol) == 0

    return pa.Check(
        ok, element_wise=True, name="neutral", description="a net formal charge", **options
    )


def stereo_defined(**options: Any) -> pa.Check:
    """Every stereocentre and stereogenic double bond has its configuration given.

    Often used as a warning rather than an error: ``stereo_defined(raise_warning=True)``.

    Args:
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """

    def ok(smiles: str) -> bool:
        """Whether one SMILES specifies all its stereochemistry (unreadable SMILES pass).

        Args:
            smiles: The SMILES.

        Returns:
            True if fully specified or unreadable.
        """
        mol = parse(smiles)
        if mol is None:
            return True
        return all(
            info.specified != Chem.StereoSpecified.Unspecified
            for info in Chem.FindPotentialStereo(mol)
        )

    return pa.Check(
        ok,
        element_wise=True,
        name="stereo_defined",
        description="a stereocentre or double bond with no configuration given",
        **options,
    )


def heavy_atoms(min_value: int = 1, max_value: int | None = None, **options: Any) -> pa.Check:
    """The number of heavy (non-hydrogen) atoms of the largest fragment is within limits.

    Args:
        min_value: The fewest heavy atoms allowed.
        max_value: The most heavy atoms allowed, or None for no upper limit.
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """
    top = float("inf") if max_value is None else max_value

    def ok(smiles: str) -> bool:
        """Whether one SMILES' largest fragment is within the limits (unreadable SMILES pass).

        Args:
            smiles: The SMILES.

        Returns:
            True if within the limits or unreadable.
        """
        mol = parse(smiles)
        if mol is None:
            return True
        n = max(f.GetNumHeavyAtoms() for f in Chem.GetMolFrags(mol, asMols=True))
        return min_value <= n <= top

    return pa.Check(
        ok,
        element_wise=True,
        name="heavy_atoms",
        description=f"heavy atoms outside {min_value} to {max_value}",
        statistics={"min_value": min_value, "max_value": max_value},
        **options,
    )


def mol_weight(max_value: float, min_value: float = 0.0, **options: Any) -> pa.Check:
    """The molecular weight (g/mol) of the largest fragment is within limits.

    Args:
        max_value: The largest weight allowed.
        min_value: The smallest weight allowed.
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """

    def ok(smiles: str) -> bool:
        """Whether one SMILES' largest fragment weighs within the limits (unreadable pass).

        Args:
            smiles: The SMILES.

        Returns:
            True if within the limits or unreadable.
        """
        mol = parse(smiles)
        if mol is None:
            return True
        piece = max(Chem.GetMolFrags(mol, asMols=True), key=lambda f: f.GetNumHeavyAtoms())
        return min_value <= Descriptors.MolWt(piece) <= max_value

    return pa.Check(
        ok,
        element_wise=True,
        name="mol_weight",
        description=f"molecular weight outside {min_value} to {max_value} g/mol",
        statistics={"min_value": min_value, "max_value": max_value},
        **options,
    )


def no_substructures(patterns: Mapping[str, str], **options: Any) -> pa.Check:
    """None of the given substructures is present (e.g. reactive groups you want kept out).

    Args:
        patterns: Each pattern's name and its SMARTS, e.g. ``{"acyl halide": "C(=O)[Cl,Br,I]"}``.
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.

    Raises:
        ValueError: If a SMARTS cannot be read.
    """
    compiled = {}
    for label, smarts in patterns.items():
        query = Chem.MolFromSmarts(smarts)
        if query is None:
            raise ValueError(f"cannot read the SMARTS for {label!r}: {smarts}")
        compiled[label] = query

    def ok(smiles: str) -> bool:
        """Whether one SMILES contains none of the patterns (unreadable SMILES pass).

        Args:
            smiles: The SMILES.

        Returns:
            True if no pattern matches or unreadable.
        """
        mol = parse(smiles)
        return mol is None or not any(mol.HasSubstructMatch(q) for q in compiled.values())

    return pa.Check(
        ok,
        element_wise=True,
        name="no_substructures",
        description=f"contains one of: {', '.join(compiled)}",
        statistics={"patterns": dict(patterns)},
        **options,
    )


def unique_molecules(**options: Any) -> pa.Check:
    """No molecule appears twice, even when written as different SMILES.

    Molecules are compared by standard InChIKey, so ``OCC`` and ``CCO`` are the same. The first
    record of each molecule passes; every later one fails. Unreadable SMILES pass.

    Args:
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """

    def ok(smiles: pd.Series) -> pd.Series:
        """Which records are the first (or only) record of their molecule.

        Args:
            smiles: The column of SMILES.

        Returns:
            True for first occurrences and unreadable SMILES.
        """
        keys = smiles.map(inchikey)
        return ~(keys.notna() & keys.duplicated(keep="first"))

    return pa.Check(
        ok,
        name="unique_molecules",
        description="the same molecule as an earlier record",
        **options,
    )


def inchikey_format(**options: Any) -> pa.Check:
    """The value is shaped like a standard InChIKey (14 letters, 10 letters, 1 letter).

    Args:
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """
    return pa.Check(
        lambda value: isinstance(value, str) and INCHIKEY.match(value) is not None,
        element_wise=True,
        name="inchikey_format",
        description="not shaped like a standard InChIKey",
        **options,
    )


def inchikey_matches(smiles_column: str, inchikey_column: str, **options: Any) -> pa.Check:
    """A table-wide check: each record's InChIKey is the one its SMILES gives.

    Records whose SMILES cannot be read, or whose InChIKey is malformed, pass: the column
    checks report those.

    Args:
        smiles_column: The name of the SMILES column.
        inchikey_column: The name of the InChIKey column.
        **options: Passed on to ``pandera.Check``.

    Returns:
        The check, to put in ``DataFrameSchema(checks=[...])``.
    """

    def ok(df: pd.DataFrame) -> pd.Series:
        """Which records' InChIKeys agree with their SMILES.

        Args:
            df: The table.

        Returns:
            True where they agree, or where either value cannot be compared.
        """
        expected = df[smiles_column].map(inchikey)
        given = df[inchikey_column]
        comparable = expected.notna() & given.astype(str).str.match(INCHIKEY.pattern)
        return ~comparable | (expected == given)

    return pa.Check(
        ok,
        name="inchikey_matches",
        description=f"{inchikey_column} is not the InChIKey of {smiles_column}",
        **options,
    )
