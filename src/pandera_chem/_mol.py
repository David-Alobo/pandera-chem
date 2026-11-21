"""Reading molecules once: a cached parser shared by every check.

A schema often runs several chemistry checks on the same column, and a table often repeats the
same SMILES. Parsing is the slow part, so each distinct SMILES is parsed and sanitised once and
the result is reused. The cached molecules are only read, never changed, by the checks.
"""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem, rdBase

CACHE_SIZE = 65_536  # distinct SMILES kept; about 100-200 MB at drug-like sizes


@lru_cache(maxsize=CACHE_SIZE)
def parse(smiles: str) -> Chem.Mol | None:
    """Read and sanitise one SMILES, quietly.

    Args:
        smiles: The SMILES to read.

    Returns:
        The sanitised molecule, or None if it cannot be read or sanitised. RDKit's own error
        messages are not printed; use :func:`pandera_chem.problems` to see why one failed.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    with rdBase.BlockLogs():
        return Chem.MolFromSmiles(smiles)


@lru_cache(maxsize=CACHE_SIZE)
def inchikey(smiles: str) -> str | None:
    """The standard InChIKey of a SMILES, or None if it cannot be read or has no key.

    Args:
        smiles: The SMILES.

    Returns:
        The InChIKey, or None.
    """
    mol = parse(smiles)
    if mol is None:
        return None
    with rdBase.BlockLogs():
        key = Chem.MolToInchiKey(mol)
    return key or None


def problems(smiles: str) -> list[str]:
    """Why a SMILES cannot be read or sanitised, in RDKit's words.

    Args:
        smiles: The SMILES.

    Returns:
        One message per problem; empty if the SMILES reads and sanitises.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return ["empty SMILES"]
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return ["the SMILES could not be read"]
        return [p.Message() for p in Chem.DetectChemistryProblems(mol)]


def clear_cache() -> None:
    """Forget every parsed molecule and InChIKey (to free memory after a large table).

    Returns:
        None.
    """
    parse.cache_clear()
    inchikey.cache_clear()
