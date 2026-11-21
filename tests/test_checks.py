import warnings

import pandas as pd
import pandera.pandas as pa
import pytest
from pandera.errors import SchemaError, SchemaWarning

import pandera_chem as pc

UNREADABLE = "c1ccccc"  # a ring that never closes


def passes(check: pa.Check, values: list) -> list[bool]:
    """Which values pass one column check, as a list of booleans."""
    schema = pa.DataFrameSchema({"x": pa.Column(object, check, nullable=True)})
    df = pd.DataFrame({"x": values})
    try:
        schema.validate(df, lazy=True)
        return [True] * len(values)
    except pa.errors.SchemaErrors as error:
        failed = set(error.failure_cases["index"])
        return [i not in failed for i in range(len(values))]


@pytest.mark.parametrize(
    ("smiles", "ok"),
    [
        ("CCO", True),
        ("c1ccccc1", True),
        (UNREADABLE, False),
        ("C(C)(C)(C)(C)C", False),  # a carbon with five bonds
        ("c1cccc1", False),  # an aromatic ring that cannot be kekulised
        ("", False),
        ("not a smiles", False),
    ],
)
def test_valid_smiles(smiles, ok):
    assert passes(pc.valid_smiles(), [smiles]) == [ok]


def test_allowed_elements_default_and_custom():
    values = ["CCO", "C[Se]C", "CC(=O)[O-].[Na+]", UNREADABLE]
    assert passes(pc.allowed_elements(), values) == [True, False, False, True]
    assert passes(pc.allowed_elements({"C", "O", "Se"}), values) == [True, True, False, True]


def test_single_fragment():
    values = ["CCO", "CC(=O)[O-].[Na+]", "c1ccccc1.O", UNREADABLE]
    assert passes(pc.single_fragment(), values) == [True, False, False, True]


def test_neutral_allows_charges_that_cancel():
    values = ["C[N+](=O)[O-]", "[NH3+]CC(=O)[O-]", "CC[NH3+]", "CC(=O)[O-]", UNREADABLE]
    assert passes(pc.neutral(), values) == [True, True, False, False, True]


def test_stereo_defined():
    values = ["C[C@H](N)C(=O)O", "CC(N)C(=O)O", "C/C=C/C", "CC=CC", "CCO", UNREADABLE]
    assert passes(pc.stereo_defined(), values) == [True, False, True, False, True, True]


def test_heavy_atoms_counts_the_largest_fragment():
    values = ["CCO", "CCCCCC", "CCO.[Na+]", "C" * 60, UNREADABLE]
    assert passes(pc.heavy_atoms(3, 50), values) == [True, True, True, False, True]
    assert passes(pc.heavy_atoms(4), values) == [False, True, False, True, True]


def test_mol_weight():
    values = ["CCO", "CC(=O)Oc1ccccc1C(=O)O", UNREADABLE]  # 46.07 and 180.16 g/mol
    assert passes(pc.mol_weight(100), values) == [True, False, True]
    assert passes(pc.mol_weight(200, min_value=50), values) == [False, True, True]


def test_no_substructures():
    check = pc.no_substructures({"acyl halide": "C(=O)[Cl,Br,I]", "aldehyde": "[CX3H1](=O)[#6]"})
    values = ["CC(=O)Cl", "CC=O", "CC(=O)O", UNREADABLE]
    assert passes(check, values) == [False, False, True, True]


def test_no_substructures_rejects_a_bad_smarts():
    with pytest.raises(ValueError, match="cannot read the SMARTS"):
        pc.no_substructures({"broken": "C(("})


def test_unique_molecules_compares_molecules_not_strings():
    values = ["CCO", "OCC", "c1ccccc1", "C1=CC=CC=C1", "CCN", UNREADABLE, UNREADABLE]
    assert passes(pc.unique_molecules(), values) == [True, False, True, False, True, True, True]


def test_inchikey_format():
    values = [
        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "LFQSCWFLJHTTHZ-UHFFFAOYSA",
        "lfqscwfljhtthz-uhfffaoysa-n",
    ]
    assert passes(pc.inchikey_format(), values) == [True, False, False]


def test_inchikey_matches():
    df = pd.DataFrame(
        {
            "smiles": ["CCO", "CCO", UNREADABLE, "CCN"],
            "key": [
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",  # ethanol: right
                "QUSNBJAOOMFDIB-UHFFFAOYSA-N",  # ethylamine's key on ethanol: wrong
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",  # SMILES unreadable: not compared
                "not a key",  # malformed: not compared
            ],
        }
    )
    schema = pa.DataFrameSchema(checks=[pc.inchikey_matches("smiles", "key")])
    with pytest.raises(pa.errors.SchemaErrors) as error:
        schema.validate(df, lazy=True)
    # a table-wide check is reported once per column; the row is the same
    assert set(error.value.failure_cases["index"]) == {1}


def test_checks_are_named_and_described():
    check = pc.valid_smiles()
    assert check.name == "valid_smiles"
    assert "RDKit" in check.description


def test_options_are_passed_to_pandera():
    schema = pa.DataFrameSchema({"x": pa.Column(str, pc.stereo_defined(raise_warning=True))})
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        schema.validate(pd.DataFrame({"x": ["CC(N)C(=O)O"]}))
    assert any(issubclass(w.category, SchemaWarning) for w in caught)


def test_a_failing_check_raises_in_plain_pandera():
    schema = pa.DataFrameSchema({"x": pa.Column(str, pc.valid_smiles())})
    with pytest.raises(SchemaError, match="valid_smiles"):
        schema.validate(pd.DataFrame({"x": ["CCO", UNREADABLE]}))


def test_problems_explains_failures():
    assert pc.problems("CCO") == []
    assert pc.problems(UNREADABLE) == ["the SMILES could not be read"]
    assert "valence" in pc.problems("C(C)(C)(C)(C)C")[0].lower()
    assert pc.problems("") == ["empty SMILES"]


def test_inchikey_helper_and_cache():
    assert pc.inchikey("CCO") == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert pc.inchikey(UNREADABLE) is None
    pc.clear_cache()
    assert pc.inchikey("OCC") == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
