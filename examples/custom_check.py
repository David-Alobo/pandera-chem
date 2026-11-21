"""Usage guide: a check of your own, written the same way as pandera-chem's."""

import pandas as pd
import pandera.pandas as pa
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors

import pandera_chem as pc


def max_rotatable_bonds(limit: int) -> pa.Check:
    """At most `limit` rotatable bonds; unreadable SMILES pass (valid_smiles reports them)."""

    def ok(smiles: str) -> bool:
        with rdBase.BlockLogs():  # keep RDKit quiet about SMILES it cannot read
            mol = Chem.MolFromSmiles(smiles)
        return mol is None or Descriptors.NumRotatableBonds(mol) <= limit

    return pa.Check(ok, element_wise=True, name="max_rotatable_bonds")


schema = pa.DataFrameSchema({"smiles": pa.Column(str, [pc.valid_smiles(), max_rotatable_bonds(5)])})
df = pd.DataFrame({"smiles": ["CCO", "CCCCCCCCCCCC", "c1ccccc"]})
print(pc.split(df, schema).quarantined.to_string())
