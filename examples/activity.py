"""Usage guide: checking activity data (values, units and pChEMBL)."""

import pandas as pd
import pandera.pandas as pa

import pandera_chem as pc

measurements = pd.DataFrame(
    {
        "smiles": ["CC(=O)Nc1ccc(O)cc1"] * 7,
        "value": [50, 0.05, 4.2, 50, 1e-6, 10, 250],
        "units": ["nM", "µM", "uM", "mM", "nM", "nm", "ug/mL"],
        "relation": ["=", "=", ">", "=", "=", "=", "="],
        "pchembl": [7.30, 7.30, None, 7.30, None, 8.00, 6.00],
    }
)
schema = pa.DataFrameSchema(
    {
        "smiles": pa.Column(str, pc.valid_smiles()),
        "value": pa.Column(float),
        "units": pa.Column(str, pc.concentration_units()),
        "relation": pa.Column(str, pa.Check.isin(["=", "<", ">"])),
        "pchembl": pa.Column(float, nullable=True),
    },
    checks=[
        pc.pchembl_matches("value", "units", "pchembl", relation_column="relation"),
        pc.activity_in_range("value", "units"),  # 1 pM to 1 M once in nM
    ],
)

result = pc.split(measurements, schema)
print(result.quarantined[["value", "units", "pchembl", "reasons"]].to_string())
