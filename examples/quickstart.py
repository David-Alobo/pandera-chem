"""The README's example: check a small table of measured compounds, and quarantine the bad rows."""

import pandas as pd
import pandera.pandas as pa

import pandera_chem as pc

table = pd.DataFrame(
    {
        "compound_id": ["C1", "C2", "C3", "C4", "C5", "C6", "C7"],
        "smiles": [
            "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
            "OC(=O)c1ccccc1OC(C)=O",  # aspirin again, written differently
            "c1ccc2ccccc2c1",  # naphthalene
            "c1ccccc",  # a ring that never closes
            "C[Se]c1ccccc1",  # selenium
            "CC(=O)[O-].[Na+]",  # a salt
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # ibuprofen, stereocentre not given
        ],
        "pic50": [5.2, 5.4, 4.1, 6.0, 5.0, 3.2, 15.0],
    }
)

schema = pa.DataFrameSchema(
    {
        "compound_id": pa.Column(str, unique=True),
        "smiles": pa.Column(
            str,
            [
                pc.valid_smiles(),
                pc.allowed_elements(),
                pc.single_fragment(),
                pc.unique_molecules(),
                pc.stereo_defined(raise_warning=True),  # a warning, not a failure
            ],
        ),
        "pic50": pa.Column(float, pa.Check.in_range(0, 14)),
    }
)

result = pc.split(table, schema)
print(result.report)
print()
print(result.quarantined[["compound_id", "reasons"]].to_string(index=False))
print()
print(pc.problems("c1ccccc"))
