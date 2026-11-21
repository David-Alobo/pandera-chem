"""Usage guide: a table that stores InChIKeys next to its SMILES."""

import pandas as pd
import pandera.pandas as pa

import pandera_chem as pc

schema = pa.DataFrameSchema(
    {
        "smiles": pa.Column(str, pc.valid_smiles()),
        "inchikey": pa.Column(str, pc.inchikey_format()),
    },
    checks=[pc.inchikey_matches("smiles", "inchikey")],  # a table-wide check
)

table = pd.DataFrame(
    {
        "smiles": ["CCO", "CCN", "c1ccccc1"],
        "inchikey": [
            pc.inchikey("CCO"),  # right
            pc.inchikey("CCO"),  # ethanol's key on ethylamine: a copy-paste slip
            "UHOVQNZJYSORNB-UHFFFAOYSA",  # benzene's key, cut short
        ],
    }
)
result = pc.split(table, schema)
print(result.quarantined["reasons"].to_string())
