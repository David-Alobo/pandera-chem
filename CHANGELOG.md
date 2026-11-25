# Changelog

## 0.4.0 (unreleased)

### Added

- **Pictures in the report card.** Each quarantined row is now a card: the molecule drawn with
  RDKit (inline SVG, still one self-contained file), the atoms or bonds that failed a check
  highlighted (the disallowed element, the extra fragment, the charged atoms, the matched
  forbidden group, the stereocentres and double bonds with no configuration), a short caption
  ("element not allowed: Se", "same molecule as row 0", or RDKit's reason for an unreadable
  SMILES), the reasons as coloured labels and the row's values. Two cards per line on a laptop,
  one on a phone. Turn off with `report_html(result, pictures=False)` or `--no-pictures`.
- `Split.schema`: the schema `split()` used, so the report can read each check's settings.

## 0.3.0 (25 November 2026)

### Added

- **Activity-data checks**, with the optional extra `pip install "pandera-chem[units]"` (pint):
  - `concentration_units()`: the unit is a molar concentration, however it is spelled (nM, uM,
    µM with either mu, mol/L...); fails ug/mL, nm, %, n/a;
  - `pchembl_matches(value, units, pchembl, relation_column=None, tolerance=0.01)`: table-wide,
    the stored pChEMBL equals 9 - log10(value in nM); with a relation column, only "=" rows;
  - `activity_in_range(value, units, min_nm=0.001, max_nm=1e9)`: table-wide, the value in nM is
    plausible and above zero (catches unit slips).
  All three can be named in rules files. Without pint they say how to install it; the rest of the
  package works as before.
- Usage guide section 10, with a tested example.

## 0.2.0 (21 November 2026)

### Added

- **Rules files**: `load_schema("rules.yaml")` and `schema_from_dict()` build a Pandera schema from
  YAML, naming pandera-chem's checks and a set of Pandera's own. A mistake (an unknown check, an
  argument a check does not take, a wrong dtype) stops with a message that says where.
- **Report cards**: `report_html()` writes one self-contained web page (counts, reasons with bars,
  warnings, the quarantined rows; light and dark); `report_json()` and `report_dict()` give the same
  for programs.
- **Command line**: `pandera-chem check data.csv --rules rules.yaml` writes the valid rows, the
  quarantined rows and the report card, and exits 0, 1 (more quarantined than `--max-quarantined`)
  or 2 (unusable input), so it can gate CI. Reads CSV, TSV, Parquet and SDF (no SDF record is
  dropped, and each keeps its data fields). `pandera-chem init` writes a starter rules file.
- Usage guide sections on rules files, the command line and a GitHub Actions job; every file and
  command shown is checked by the tests.

### Changed

- Depends on PyYAML.

## 0.1.1 (21 November 2026)

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
