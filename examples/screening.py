"""Usage guide: rules for a screening library."""

import pandas as pd
import pandera.pandas as pa

import pandera_chem as pc

REACTIVE = {
    "acyl halide": "C(=O)[Cl,Br,I]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "isocyanate": "N=C=O",
}
library = pa.DataFrameSchema(
    {
        "smiles": pa.Column(
            str,
            [
                pc.valid_smiles(),
                pc.allowed_elements({"H", "C", "N", "O", "F", "S", "Cl", "Br"}),
                pc.single_fragment(),
                pc.neutral(),
                pc.heavy_atoms(6, 40),
                pc.mol_weight(500),
                pc.no_substructures(REACTIVE),
                pc.unique_molecules(),
            ],
        )
    }
)

candidates = pd.DataFrame(
    {
        "smiles": [
            "CC(=O)Nc1ccc(O)cc1",  # paracetamol
            "O=Cc1ccc(N)cc1",  # an aldehyde
            "CCO",  # too small
            "CC(=O)Nc1ccc(OP(=O)(O)O)cc1",  # phosphorus: not in this library's list
            "Cc1ccccc1N=C=O",  # an isocyanate
            "CN1CCC[C@H]1c1cccnc1",  # nicotine
        ]
    }
)
result = pc.split(candidates, library)
print(result.quarantined.to_string())
