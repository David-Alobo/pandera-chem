# pandera-chem

[![tests](https://github.com/David-Alobo/pandera-chem/actions/workflows/ci.yml/badge.svg)](https://github.com/David-Alobo/pandera-chem/actions/workflows/ci.yml)
![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776ab)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

**Chemistry checks for [Pandera](https://pandera.readthedocs.io/).** Validate tables of molecules
with [RDKit](https://www.rdkit.org/) inside the data pipeline you already have, and set the bad rows
aside with the reason each one failed, instead of crashing or losing them.

Pandera checks the columns of a table: types, ranges, uniqueness. pandera-chem adds checks that
understand what is *in* a SMILES column: whether each molecule can be read, which elements it
contains, whether it is a salt or a mixture, whether it is charged, whether its stereochemistry is
given, and whether the same molecule appears twice under different SMILES.

## Install

```bash
pip install pandera-chem
```

or, with [uv](https://docs.astral.sh/uv/): `uv add pandera-chem`. It needs Python 3.11 or newer, and
installs Pandera and RDKit for you.

## Example

A table of measured compounds, with a few of the problems real data has:

```python
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
```

It prints:

```text
7 rows in: 2 valid, 5 quarantined
       2  smiles: allowed_elements
       1  smiles: valid_smiles
       1  smiles: single_fragment
       1  smiles: unique_molecules
       1  pic50: in_range(0, 14)
  warning: Column 'smiles' failed element-wise validator number 4: <Check stereo_defined> failure cases: CC(C)Cc1ccc(cc1)C(C)C(=O)O

compound_id                                           reasons
         C2                          smiles: unique_molecules
         C4                              smiles: valid_smiles
         C5                          smiles: allowed_elements
         C6 smiles: allowed_elements; smiles: single_fragment
         C7                            pic50: in_range(0, 14)

['the SMILES could not be read']
```

- `result.valid` holds the rows that passed (aspirin and naphthalene), ready for the next step.
- `result.quarantined` holds the rest, unchanged, with a `reasons` column: save it, fix the records
  and run them again.
- `result.report` has the counts, and `result.report.reconciles` is `True` when every row is
  accounted for (rows in = valid + quarantined). Log it with every run.
- Ibuprofen passed with a warning: its stereocentre is not given, which you may accept for a
  racemate. Any check can be made a warning with `raise_warning=True`.

The checks also work in plain Pandera: `schema.validate(table)` raises an error naming the check that
failed, as Pandera always does.

## The checks

Put these on a column of SMILES (the last two are for InChIKeys):

| Check | Fails when |
|---|---|
| `valid_smiles()` | RDKit cannot read or sanitise the SMILES (syntax, valence, aromaticity) |
| `allowed_elements(elements=ORGANIC)` | an atom is not in the allowed set (default: H, B, C, N, O, F, P, S, Cl, Br, I) |
| `single_fragment()` | the record holds more than one molecule: a salt, a solvent or a mixture |
| `neutral()` | the net formal charge is not zero (charges that cancel, as in a nitro group, are fine) |
| `stereo_defined()` | a stereocentre or stereogenic double bond has no configuration given |
| `heavy_atoms(min_value, max_value)` | the largest fragment has too few or too many heavy atoms |
| `mol_weight(max_value, min_value=0)` | the largest fragment is too light or too heavy (g/mol) |
| `no_substructures({"name": "SMARTS", ...})` | the molecule contains one of the given groups |
| `unique_molecules()` | the same molecule appeared in an earlier row, even as a different SMILES |
| `inchikey_format()` | (InChIKey column) the value is not shaped like a standard InChIKey |
| `inchikey_matches("smiles", "inchikey")` | (table-wide) a row's InChIKey is not the one its SMILES gives |

**One reason per problem.** Only `valid_smiles()` fails a SMILES that cannot be read; every other
check lets it pass. So a broken record is reported once, as unreadable, not under every check.

Every check accepts Pandera's own options, for example `raise_warning=True` or `n_failure_cases=10`.
Each distinct SMILES is parsed once and reused by all the checks; call `pc.clear_cache()` after a very
large table to free the memory.

`pc.problems(smiles)` says why a SMILES cannot be read, in RDKit's words, and `pc.inchikey(smiles)`
gives its standard InChIKey.

## What it is not

pandera-chem checks molecules; it does not change them. To repair records (strip salts, neutralise
charges, pick one tautomer), use RDKit's
[MolStandardize](https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html)
or the [ChEMBL structure pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline) first, then
check the result with pandera-chem. For a web page that validates an uploaded file, see
[ChemAudit](https://nfdi4chem.de/introducing-chemaudit/).

## Learn the ideas behind it

pandera-chem grew out of the free course
[Data Engineering for Cheminformatics](https://monstrous.com.ng/learn/), which builds these checks
by hand in Module 6 (validating molecules) and Module 12 (validating data with Pydantic and
Pandera), and the quarantine in Module 10 (pipelines that survive failure).

## Contributing

Issues and pull requests are welcome. To work on it:

```bash
git clone https://github.com/David-Alobo/pandera-chem.git
cd pandera-chem
uv sync
uv run pytest
uv run ruff check . && uv run ruff format --check .
```

The example above is run by the tests, which check that this README shows what it really prints.

## Licence

MIT ([LICENSE](LICENSE)). Made by [Monstrous Data](https://monstrous.com.ng).
