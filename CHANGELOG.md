# Changelog

## 0.1.0 (unreleased)

First version.

- Checks for a column of SMILES: `valid_smiles`, `allowed_elements`, `single_fragment`, `neutral`,
  `stereo_defined`, `heavy_atoms`, `mol_weight`, `no_substructures`, `unique_molecules`.
- Checks for InChIKeys: `inchikey_format`, and the table-wide `inchikey_matches`.
- `split()`: the rows that pass, the rows that fail with their reasons, and a report whose counts
  reconcile.
- `problems()` and `inchikey()` helpers; each distinct SMILES is parsed once and cached.
