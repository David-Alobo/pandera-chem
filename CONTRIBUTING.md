# Contributing to pandera-chem

Thank you for helping. The most useful contributions are:

- **A wrong answer from a check.** A SMILES that a check passes or fails when it should not. Open a
  [bug report](https://github.com/David-Alobo/pandera-chem/issues/new?template=bug.yml) with the
  SMILES and what you expected.
- **A check you need.** Describe the problem in your data in a
  [suggestion](https://github.com/David-Alobo/pandera-chem/issues/new?template=feature.yml), with
  examples that should pass and fail.
- **A pull request.** See below.

## Setting up

You need [uv](https://docs.astral.sh/uv/) and Git.

```bash
git clone https://github.com/David-Alobo/pandera-chem.git
cd pandera-chem
uv sync                      # makes .venv with everything, from uv.lock
uv run pytest --cov          # the tests, with coverage
uv run ruff check . && uv run ruff format --check .
```

Optionally, run the same checks before every commit:
`uv tool install pre-commit && pre-commit install`.

## How the code is organised

- `src/pandera_chem/checks.py`: the checks. Each is a function that returns a `pandera.Check`.
- `src/pandera_chem/quarantine.py`: `split()` and its report.
- `src/pandera_chem/_mol.py`: parsing, shared and cached by every check.
- `examples/`: every example in the README and `docs/usage.md`; `tests/test_examples.py` runs each
  one and checks the page shows exactly its code and what it prints. Add a script there for any new
  example, mark it in the page with `<!-- example: examples/name.py -->`, and paste its output.

## Rules for a new check

1. **One reason per problem.** Only `valid_smiles()` fails a SMILES that cannot be read; every other
   check returns True for it, so a broken record is reported once.
2. **Parse through `_mol.parse`**, never `Chem.MolFromSmiles` directly, so each SMILES is read once.
3. **Never change the molecule.** The cached molecules are shared between checks.
4. **Accept `**options`** and pass them to `pa.Check`, and give the check a `name` (used as the reason)
   and a `description`.
5. **Test it** with molecules that pass, molecules that fail, and an unreadable SMILES.
6. Give it a Google-style docstring, add it to `__init__.py`, the README's table and `CHANGELOG.md`.

Keep each pull request to one change. By contributing you agree that your code is released under the
MIT licence. Everyone taking part follows the [code of conduct](CODE_OF_CONDUCT.md).
