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
- [7. Rules in a file](#7-rules-in-a-file)
- [8. From the command line](#8-from-the-command-line)
- [9. In continuous integration](#9-in-continuous-integration)
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

## 7. Rules in a file

Rules written in Python live in your code. Rules in a YAML file can sit next to the data, be
reviewed like any other change, and be used by people who do not write Python. A rules file has:

- `columns`: each column's name, and under it any of
  - `dtype`: `str`, `float`, `int` or `bool`;
  - `checks`: a list of checks;
  - `nullable` (default `false`), `unique` (default `false`), `required` (default `true`),
    `coerce` (default `false`) and `description`;
- `checks` (optional): table-wide checks, such as `inchikey_matches`;
- `strict` (optional): `true` to refuse columns the rules do not name.

Each check is its **name** (`valid_smiles`) or its **name with its arguments**
(`heavy_atoms: {min_value: 6, max_value: 50}`), using the same arguments as in Python, including
Pandera's options such as `raise_warning: true`. Pandera's own `in_range`, `isin`, `notin`, `gt`, `ge`,
`lt`, `le`, `eq`, `ne` and `str_matches`/`str_contains`/`str_startswith`/`str_endswith`/`str_length`
can be used too. Here is the rules file used in the [README](../README.md#rules-files-and-the-command-line):

<!-- file: examples/rules.yaml -->
```yaml
# Rules for a table of measured compounds.
columns:
  compound_id:
    dtype: str
    unique: true
  smiles:
    dtype: str
    checks:
      - valid_smiles
      - allowed_elements
      - single_fragment
      - unique_molecules
      - stereo_defined: {raise_warning: true}
  pic50:
    dtype: float
    checks:
      - in_range: {min_value: 0, max_value: 14}
```

`pc.load_schema()` turns it into an ordinary Pandera schema, and `pc.report_html()` and
`pc.report_json()` write the report card:

<!-- example: examples/rules_file.py -->
```python
from pathlib import Path

import pandas as pd

import pandera_chem as pc

HERE = Path(__file__).parent
schema = pc.load_schema(HERE / "rules.yaml")
table = pd.read_csv(HERE / "compounds.csv")

result = pc.split(table, schema)
Path("report.html").write_text(pc.report_html(result, source="compounds.csv"), encoding="utf-8")
Path("report.json").write_text(pc.report_json(result, source="compounds.csv"), encoding="utf-8")

print(result.report.rows_valid, "valid rows;", result.report.rows_quarantined, "quarantined")
print(pc.report_dict(result)["by_reason"])
```

```text
2 valid rows; 5 quarantined
{'smiles: allowed_elements': 2, 'smiles: valid_smiles': 1, 'smiles: single_fragment': 1, 'smiles: unique_molecules': 1, 'pic50: in_range(0, 14)': 1}
```

A mistake in the rules stops before any data is checked, with a message that says where: a check
name that does not exist (`unknown check 'valid_smile'`), an argument it does not take
(`heavy_atoms: unknown arguments ['min']`), or a dtype that is not one of the four.

## 8. From the command line

Installing pandera-chem also installs a `pandera-chem` command (or run `python -m pandera_chem`).
Start with a rules file to edit:

```bash
pandera-chem init                # writes rules.yaml, with a comment on every line
```

Then check a data file against it:

<!-- command: pandera-chem check examples/compounds.csv --rules examples/rules.yaml --out report -->
```bash
pandera-chem check examples/compounds.csv --rules examples/rules.yaml --out report
```

```text
7 rows in: 2 valid, 5 quarantined
       2  smiles: allowed_elements
       1  smiles: valid_smiles
       1  smiles: single_fragment
       1  smiles: unique_molecules
       1  pic50: in_range(0, 14)
  warning: Column 'smiles' failed element-wise validator number 4: <Check stereo_defined> failure cases: CC(C)Cc1ccc(cc1)C(C)C(=O)O
wrote report/compounds.{valid.csv,quarantine.csv,report.html,report.json}
FAIL: 71.4% of rows quarantined (allowed: 0.0%)
```

In the output folder (`--out`, default `pandera-chem-out`):

| File | What it holds |
|---|---|
| `compounds.valid.csv` | the rows that passed |
| `compounds.quarantine.csv` | the rows that failed, with a `reasons` column |
| `compounds.report.html` | the report card: open it in a browser, attach it to a CI run or email it |
| `compounds.report.json` | the same counts for programs |

**Exit codes**, for scripts and CI: `0` when the share of quarantined rows is within
`--max-quarantined`, `1` when it is above it, and `2` when the file or the rules cannot be used
(a missing file, a mistake in the rules, a column the rules need but the file lacks). The default
`--max-quarantined` is `0`: any bad row fails the run. `--max-quarantined 0.05` allows up to 5%.

**Files it reads:**

- `.csv`, and `.tsv` or `.txt` (tab-separated). Every column is read as text, so identifiers keep
  their leading zeros; columns the rules give a `float` or `int` dtype are converted to numbers,
  and text that is not a number (such as `n/a`) becomes empty, so the column's `nullable` setting
  decides whether that row is quarantined.
- `.parquet`, with `pyarrow` installed (`pip install pyarrow`).
- `.sdf`: one row per record, its SMILES in a `smiles` column (or `--smiles-column`) and its data
  fields as the other columns. No record is dropped: one RDKit cannot read keeps its data fields and
  an empty or unsanitised SMILES, so `valid_smiles` quarantines it and you can see which it was.

## 9. In continuous integration

Because the command fails with exit code 1, it can guard a repository of data: a pull request that
adds bad records fails, and the report card is kept with the run. A GitHub Actions job:

```yaml
# .github/workflows/data.yml
name: data quality
on: [push, pull_request]
jobs:
  check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v7
      - uses: actions/setup-python@v6
        with:
          python-version: "3.12"
      - run: pip install pandera-chem
      - run: pandera-chem check data/compounds.csv --rules data/rules.yaml --out report
      - uses: actions/upload-artifact@v4
        if: always()   # keep the report card when the check fails, too
        with:
          name: data-quality-report
          path: report/
```

To allow a small share of bad rows while they are being fixed, add `--max-quarantined 0.02`.

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
