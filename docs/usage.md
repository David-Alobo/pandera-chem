# Using pandera-chem

A practical guide, from installing it to running it in a pipeline. The [README](../README.md) has
the first example and the list of checks; this guide shows how to use them for real work. Every
example here is a script in [`examples/`](../examples/), and the tests run each one and check that
this page shows exactly the code and what it prints.

- [Install](#install)
- [1. Stop when anything is wrong](#1-stop-when-anything-is-wrong)
- [2. Keep going: quarantine the bad rows](#2-keep-going-quarantine-the-bad-rows)
- [3. In a pipeline](#3-in-a-pipeline)
- [4. Rules for a screening library](#4-rules-for-a-screening-library)
- [5. A table that stores InChIKeys](#5-a-table-that-stores-inchikeys)
- [6. A check of your own](#6-a-check-of-your-own)
- [Things to know](#things-to-know)

## Install

```bash
pip install pandera-chem
```

or `uv add pandera-chem`. In a Jupyter or Google Colab notebook, run this in a cell:

```text
!pip install pandera-chem
```

It needs Python 3.11 or newer, and installs Pandera and RDKit.

## 1. Stop when anything is wrong

A Pandera *schema* lists each column and the checks it must pass; pandera-chem's checks go in the
same list as Pandera's own. `schema.validate()` returns the table if every row passes, and raises an
error if any fails. With `lazy=True` the error lists every failure, not just the first, and its
`failure_cases` table says which row, column and check:

<!-- example: examples/strict.py -->
```python
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
```

```text
 index column            check failure_case
     1 smiles     valid_smiles      c1ccccc
     3 smiles allowed_elements       C[Se]C
     2 smiles  single_fragment       CCN.Cl
```

Use this when bad data should stop everything: in a test, or when a file is checked before it is
accepted from a supplier.

## 2. Keep going: quarantine the bad rows

`pc.split(df, schema)` runs the same schema but does not stop. It returns three things:

| | What it holds |
|---|---|
| `result.valid` | the rows that passed, ready for the next step |
| `result.quarantined` | the rows that failed, unchanged, with a `reasons` column |
| `result.report` | the counts: rows in, valid, quarantined, and how many failed each check |

`result.report.reconciles` is `True` when every row is accounted for (rows in = valid +
quarantined). The [README's example](../README.md#example) shows it on a table of seven compounds.

## 3. In a pipeline

Call `split()` where data enters your pipeline, after reading it and before anything else uses it.
Then pass the valid rows on, save the quarantined ones where someone will fix them, log the report,
and stop the run if far more rows fail than usual, which usually means the file itself is wrong:

<!-- example: examples/pipeline.py -->
```python
import datetime
import logging
import sys
from pathlib import Path

import pandas as pd
import pandera.pandas as pa

import pandera_chem as pc

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s", stream=sys.stdout)
log = logging.getLogger("ingest")

SCHEMA = pa.DataFrameSchema(
    {
        "compound_id": pa.Column(str, unique=True),
        "smiles": pa.Column(str, [pc.valid_smiles(), pc.single_fragment(), pc.unique_molecules()]),
        "pic50": pa.Column(float, pa.Check.in_range(0, 14)),
    }
)
MAX_QUARANTINED = 0.25  # a larger share usually means something is wrong upstream


def ingest(path: Path, out: Path) -> pd.DataFrame:
    """Read a supplier's file, keep the rows that pass, and save the others with their reasons."""
    df = pd.read_csv(path).reset_index(drop=True)  # split() needs a unique index
    result = pc.split(df, SCHEMA)
    log.info("%s\n%s", path.name, result.report)
    if not result.report.reconciles:
        raise RuntimeError("rows went missing")
    stamp = datetime.date(2026, 11, 21).isoformat()  # in a real job: datetime.date.today()
    out.mkdir(exist_ok=True)
    result.valid.to_csv(out / f"valid-{stamp}.csv", index=False)
    result.quarantined.to_csv(out / f"quarantine-{stamp}.csv", index=False)
    share = result.report.rows_quarantined / result.report.rows_in
    if share > MAX_QUARANTINED:
        raise RuntimeError(f"{share:.0%} of rows quarantined; check the supplier's file")
    return result.valid


# a supplier's file, as it might arrive
pd.DataFrame(
    {
        "compound_id": ["S-1", "S-2", "S-3", "S-4", "S-5", "S-6", "S-7", "S-8"],
        "smiles": ["CCO", "CCN", "c1ccccc1O", "CCN.Cl", "CC(=O)O", "OCC", "CCCC", "CCOC"],
        "pic50": [5.1, 6.3, 4.8, 5.5, 7.0, 5.1, 4.2, 6.6],
    }
).to_csv("supplier.csv", index=False)

clean = ingest(Path("supplier.csv"), Path("out"))
print(len(clean), "rows go on to the next step")
print(sorted(p.name for p in Path("out").iterdir()))
```

```text
INFO supplier.csv
8 rows in: 6 valid, 2 quarantined
       1  smiles: single_fragment
       1  smiles: unique_molecules
6 rows go on to the next step
['quarantine-2026-11-21.csv', 'valid-2026-11-21.csv']
```

`OCC` was set aside as a second copy of `CCO` (both are ethanol), and `CCN.Cl` as a salt. The two
CSV files hold the rows and their reasons; the quarantine file can be fixed and fed back in.

## 4. Rules for a screening library

A library for screening usually wants drug-like sizes, a short list of elements, and no reactive
groups. `no_substructures()` takes the groups as SMARTS patterns, each with a name:

<!-- example: examples/screening.py -->
```python
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
```

```text
                        smiles                   reasons
1               O=Cc1ccc(N)cc1  smiles: no_substructures
2                          CCO       smiles: heavy_atoms
3  CC(=O)Nc1ccc(OP(=O)(O)O)cc1  smiles: allowed_elements
4               Cc1ccccc1N=C=O  smiles: no_substructures
```

Paracetamol and nicotine pass. Your own element list replaces the default one
(`pc.ORGANIC`: H, B, C, N, O, F, P, S, Cl, Br, I).

## 5. A table that stores InChIKeys

When a table keeps an InChIKey next to each SMILES, two things can go wrong: a key that is not shaped
like an InChIKey, and a key that belongs to a different molecule (a copy-paste slip).
`inchikey_matches()` compares two columns, so it goes in the schema's table-wide `checks`:

<!-- example: examples/inchikeys.py -->
```python
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
```

```text
1             inchikey_matches
2    inchikey: inchikey_format
```

A table-wide check's reason has no column name. The malformed key is reported once, by
`inchikey_format()`; `inchikey_matches()` skips keys it cannot compare.

## 6. A check of your own

A check is a function that returns a `pandera.Check`. Follow the same rule as pandera-chem's checks,
letting an unreadable SMILES pass so that only `valid_smiles()` reports it, and give the check a
`name`, which `split()` uses as the reason:

<!-- example: examples/custom_check.py -->
```python
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
```

```text
         smiles                      reasons
1  CCCCCCCCCCCC  smiles: max_rotatable_bonds
2       c1ccccc         smiles: valid_smiles
```

If it would be useful to others, [suggest it](https://github.com/David-Alobo/pandera-chem/issues/new?template=feature.yml).

## Things to know

- **The index must be unique.** `split()` matches each failure to its row by the index. If yours has
  repeats, call `df.reset_index(drop=True)` first; `split()` refuses a repeated index rather than
  guess.
- **Missing columns stop everything.** A failure that belongs to no row, such as a column that is
  not there, raises Pandera's `SchemaErrors`, because the whole table is unusable.
- **Warnings do not quarantine.** A check made with `raise_warning=True` lets the row through and
  adds its message to `result.report.warnings`.
- **Empty cells.** Pandera skips empty values in a check; to forbid them, keep the column's default
  `nullable=False`, which reports them as `not_nullable`.
- **Molecules are never changed.** To repair records (strip salts, neutralise charges), use RDKit's
  MolStandardize or the ChEMBL structure pipeline first, then check the result.
- **Memory.** Each distinct SMILES is parsed once and kept for the other checks. After a very large
  table, call `pc.clear_cache()`.
- **Why a SMILES fails.** `pc.problems("C(C)(C)(C)(C)C")` gives RDKit's reason, here an impossible
  valence.
