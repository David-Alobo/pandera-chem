"""Pictures for the report card: each quarantined molecule, with what is wrong highlighted.

For every check that failed on a row's SMILES, :func:`explain` finds the atoms and bonds to
blame (the selenium atom, the counter-ion, the matched aldehyde, the stereocentre with no
configuration) and a short caption; :func:`picture` draws the molecule as an inline SVG. Nothing
is changed: the drawing works on its own copy of the molecule.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

from ._mol import inchikey, parse, problems
from .checks import ORGANIC

HIGHLIGHT = (0.97, 0.55, 0.55)  # a soft red, readable on the white tile
SIZE = (220, 160)
# the checks that look at a SMILES column (their names are the reasons in split())
MOLECULE_CHECKS = frozenset(
    {
        "valid_smiles",
        "allowed_elements",
        "single_fragment",
        "neutral",
        "stereo_defined",
        "heavy_atoms",
        "mol_weight",
        "no_substructures",
        "unique_molecules",
    }
)


@dataclass
class Explanation:
    """What to highlight in one molecule, and the words that go under its picture.

    Attributes:
        atoms: Indices of the atoms to highlight.
        bonds: Indices of the bonds to highlight.
        captions: One short sentence per failed check.
    """

    atoms: set[int] = field(default_factory=set)
    bonds: set[int] = field(default_factory=set)
    captions: list[str] = field(default_factory=list)


def _bonds_between(mol: Chem.Mol, atoms: set[int]) -> set[int]:
    """The bonds whose two atoms are both in a set.

    Args:
        mol: The molecule.
        atoms: Atom indices.

    Returns:
        The indices of the bonds inside the set.
    """
    return {b.GetIdx() for b in mol.GetBonds() if {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} <= atoms}


def _largest_fragment(mol: Chem.Mol) -> tuple[int, ...]:
    """The atom indices of the fragment with the most heavy atoms.

    Args:
        mol: The molecule.

    Returns:
        The indices of that fragment's atoms.
    """
    frags = Chem.GetMolFrags(mol)
    return max(frags, key=lambda f: sum(mol.GetAtomWithIdx(i).GetAtomicNum() > 1 for i in f))


def explain(
    smiles: str,
    failed: Iterable[str],
    statistics: Mapping[str, Mapping] | None = None,
    repeats: object = None,
) -> tuple[Chem.Mol | None, Explanation]:
    """Work out what to highlight in a molecule for the checks it failed.

    Args:
        smiles: The row's SMILES.
        failed: The names of the molecule checks the row failed (e.g. "allowed_elements").
        statistics: Each check's settings, by name, as the schema holds them (the allowed
            elements, the forbidden patterns); defaults are used where missing.
        repeats: For ``unique_molecules``, the row this one repeats, if known.

    Returns:
        The molecule (None if it cannot be read) and what to highlight and say.
    """
    statistics = statistics or {}
    why = Explanation()
    mol = parse(smiles)
    failed = list(failed)
    if mol is None:
        reason = (problems(smiles) or ["cannot be read"])[0]
        why.captions.append(f"cannot be read: {reason}")
        return None, why
    for name in failed:
        if name == "allowed_elements":
            allowed = set(statistics.get(name, {}).get("elements", ORGANIC))
            atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() not in allowed}
            symbols = sorted({mol.GetAtomWithIdx(i).GetSymbol() for i in atoms})
            why.atoms |= atoms
            why.captions.append(f"element not allowed: {', '.join(symbols)}")
        elif name == "single_fragment":
            main = set(_largest_fragment(mol))
            extra = {a.GetIdx() for a in mol.GetAtoms()} - main
            pieces = sorted(
                Chem.MolFragmentToSmiles(mol, list(f))
                for f in Chem.GetMolFrags(mol)
                if not set(f) <= main
            )
            why.atoms |= extra
            why.bonds |= _bonds_between(mol, extra)
            why.captions.append(f"extra piece: {', '.join(pieces)}")
        elif name == "neutral":
            atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge() != 0}
            why.atoms |= atoms
            why.captions.append(f"net charge {Chem.GetFormalCharge(mol):+d}")
        elif name == "stereo_defined":
            centres = bonds = 0
            for info in Chem.FindPotentialStereo(mol):
                if info.specified != Chem.StereoSpecified.Unspecified:
                    continue
                if info.type == Chem.StereoType.Atom_Tetrahedral:
                    why.atoms.add(info.centeredOn)
                    centres += 1
                elif info.type == Chem.StereoType.Bond_Double:
                    why.bonds.add(info.centeredOn)
                    bonds += 1
            parts = [f"{centres} stereocentre{'s' * (centres != 1)}"] * (centres > 0)
            parts += [f"{bonds} double bond{'s' * (bonds != 1)}"] * (bonds > 0)
            why.captions.append(f"no configuration: {' and '.join(parts) or 'stereochemistry'}")
        elif name == "no_substructures":
            patterns = statistics.get(name, {}).get("patterns", {})
            found = []
            for label, smarts in patterns.items():
                query = Chem.MolFromSmarts(smarts)
                matches = mol.GetSubstructMatches(query) if query is not None else ()
                if matches:
                    found.append(label)
                    for match in matches:
                        why.atoms |= set(match)
                        why.bonds |= _bonds_between(mol, set(match))
            why.captions.append(f"contains: {', '.join(found) or 'a forbidden group'}")
        elif name == "heavy_atoms":
            n = max(
                sum(mol.GetAtomWithIdx(i).GetAtomicNum() > 1 for i in f)
                for f in Chem.GetMolFrags(mol)
            )
            why.captions.append(f"{n} heavy atoms")
        elif name == "mol_weight":
            piece = Chem.RWMol(mol)
            keep = set(_largest_fragment(mol))
            for i in sorted({a.GetIdx() for a in mol.GetAtoms()} - keep, reverse=True):
                piece.RemoveAtom(i)
            why.captions.append(f"{Descriptors.MolWt(piece):.0f} g/mol")
        elif name == "unique_molecules":
            again = f" as row {repeats}" if repeats is not None else " as an earlier row"
            why.captions.append(f"same molecule{again}")
    return mol, why


def picture(mol: Chem.Mol, why: Explanation, size: tuple[int, int] = SIZE) -> str:
    """Draw a molecule as SVG, with the blamed atoms and bonds highlighted.

    Args:
        mol: The molecule (it is copied, not changed).
        why: What to highlight.
        size: Width and height in pixels.

    Returns:
        The SVG element, without the XML declaration, ready to put in a page.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    options = drawer.drawOptions()
    options.clearBackground = True
    options.setBackgroundColour((1, 1, 1, 1))
    atoms, bonds = sorted(why.atoms), sorted(why.bonds)
    with rdBase.BlockLogs():
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            Chem.Mol(mol),
            highlightAtoms=atoms,
            highlightBonds=bonds,
            highlightAtomColors={i: HIGHLIGHT for i in atoms},
            highlightBondColors={i: HIGHLIGHT for i in bonds},
        )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return re.sub(r"^<\?xml[^>]*\?>\s*", "", svg)


def first_rows(smiles: Mapping[object, str]) -> dict[str, object]:
    """For each molecule, the first row it appears in (by InChIKey).

    Args:
        smiles: Each row's label and SMILES, in table order.

    Returns:
        Each InChIKey and the label of its first row.
    """
    first: dict[str, object] = {}
    for label, value in smiles.items():
        key = inchikey(value) if isinstance(value, str) else None
        if key is not None and key not in first:
            first[key] = label
    return first
