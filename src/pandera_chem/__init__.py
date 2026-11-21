"""Chemistry checks for Pandera: validate tables of molecules, and quarantine the bad rows.

Quick start::

    import pandera.pandas as pa
    import pandera_chem as pc

    schema = pa.DataFrameSchema({
        "smiles": pa.Column(str, [pc.valid_smiles(), pc.allowed_elements(), pc.single_fragment()]),
        "pic50": pa.Column(float, pa.Check.in_range(0, 14)),
    })
    result = pc.split(table, schema)
    print(result.report)
"""

from ._mol import clear_cache, inchikey, problems
from .checks import (
    INCHIKEY,
    ORGANIC,
    allowed_elements,
    heavy_atoms,
    inchikey_format,
    inchikey_matches,
    mol_weight,
    neutral,
    no_substructures,
    single_fragment,
    stereo_defined,
    unique_molecules,
    valid_smiles,
)
from .quarantine import Report, Split, split

__version__ = "0.1.0"

__all__ = [
    "INCHIKEY",
    "ORGANIC",
    "Report",
    "Split",
    "allowed_elements",
    "clear_cache",
    "heavy_atoms",
    "inchikey",
    "inchikey_format",
    "inchikey_matches",
    "mol_weight",
    "neutral",
    "no_substructures",
    "problems",
    "single_fragment",
    "split",
    "stereo_defined",
    "unique_molecules",
    "valid_smiles",
]
