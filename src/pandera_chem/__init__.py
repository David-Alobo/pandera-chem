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
from ._version import __version__
from .activity import activity_in_range, concentration_units, pchembl_matches
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
from .config import RulesError, load_schema, schema_from_dict
from .quarantine import Report, Split, split
from .report import report_dict, report_html, report_json

__all__ = [
    "INCHIKEY",
    "ORGANIC",
    "Report",
    "RulesError",
    "Split",
    "activity_in_range",
    "allowed_elements",
    "clear_cache",
    "concentration_units",
    "heavy_atoms",
    "inchikey",
    "inchikey_format",
    "inchikey_matches",
    "load_schema",
    "mol_weight",
    "neutral",
    "no_substructures",
    "pchembl_matches",
    "problems",
    "report_dict",
    "report_html",
    "report_json",
    "schema_from_dict",
    "single_fragment",
    "split",
    "stereo_defined",
    "unique_molecules",
    "valid_smiles",
    "__version__",
]
