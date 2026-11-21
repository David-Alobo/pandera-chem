# Changelog

## Unreleased

### Fixed

- A table-wide check (such as `inchikey_matches`) gave one reason per column in `split()`
  ("smiles: inchikey_matches; inchikey: inchikey_matches"); it now gives one reason,
  "inchikey_matches".
- `inchikey_format()` was reported under Pandera's internal name (`str_matches(...)`); it is now
  reported as `inchikey_format`.

### Added

- A usage guide (`docs/usage.md`): checking strictly, a pipeline step that saves the quarantine and
  stops a bad run, screening-library rules, tables with InChIKeys, and writing your own check. Every
  example is a script in `examples/`, run by the tests.

### Changed

- Repository: banner, issue and pull-request forms, contributing guide, code of conduct, citation
  file, Dependabot, pre-commit, and test coverage in CI (at least 95%).

## 0.1.0 (21 November 2026)

First version.

- Checks for a column of SMILES: `valid_smiles`, `allowed_elements`, `single_fragment`, `neutral`,
  `stereo_defined`, `heavy_atoms`, `mol_weight`, `no_substructures`, `unique_molecules`.
- Checks for InChIKeys: `inchikey_format`, and the table-wide `inchikey_matches`.
- `split()`: the rows that pass, the rows that fail with their reasons, and a report whose counts
  reconcile.
- `problems()` and `inchikey()` helpers; each distinct SMILES is parsed once and cached.
