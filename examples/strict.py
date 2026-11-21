"""Usage guide: check a table strictly, and see every failure at once."""

import pandas as pd
import pandera.pandas as pa

import pandera_chem as pc

df = pd.DataFrame({"smiles": ["CCO", "c1ccccc", "CCN.Cl", "C[Se]C"]})
schema = pa.DataFrameSchema(
    {"smiles": pa.Column(str, [pc.valid_smiles(), pc.allowed_elements(), pc.single_fragment()])}
)

try:
    schema.validate(df, lazy=True)  # lazy=True: collect every failure, not just the first
except pa.errors.SchemaErrors as error:
    failures = error.failure_cases[["index", "column", "check", "failure_case"]]
    print(failures.to_string(index=False))
