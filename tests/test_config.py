import pandas as pd
import pandera.pandas as pa
import pytest

import pandera_chem as pc
from pandera_chem.config import CHEM_CHECKS, EXAMPLE_RULES

RULES = {
    "columns": {
        "compound_id": {"dtype": "str", "unique": True},
        "smiles": {
            "dtype": "str",
            "checks": [
                "valid_smiles",
                {"allowed_elements": {"elements": ["C", "H", "N", "O", "Cl"]}},
                "single_fragment",
                {"heavy_atoms": {"min_value": 2}},
                {"stereo_defined": {"raise_warning": True}},
                {"no_substructures": {"patterns": {"aldehyde": "[CX3H1](=O)[#6]"}}},
            ],
        },
        "pic50": {
            "dtype": "float",
            "nullable": True,
            "checks": [{"in_range": {"min_value": 0, "max_value": 14}}],
        },
        "inchikey": {"dtype": "str", "required": False, "checks": ["inchikey_format"]},
    },
}
TABLE = pd.DataFrame(
    {
        "compound_id": ["A", "B", "C", "D", "E"],
        "smiles": ["CCO", "c1ccccc", "CC=O", "CCN.Cl", "C[Se]C"],
        "pic50": [5.0, 6.0, None, 15.0, 4.0],
    }
)


def test_a_schema_from_rules_checks_like_one_written_in_python():
    result = pc.split(TABLE, pc.schema_from_dict(RULES))
    assert result.quarantined["reasons"].to_dict() == {
        1: "smiles: valid_smiles",
        2: "smiles: no_substructures",
        3: "smiles: single_fragment; pic50: in_range(0, 14)",
        4: "smiles: allowed_elements",
    }
    assert list(result.valid.index) == [0]


def test_a_table_wide_check_from_rules():
    schema = pc.schema_from_dict(
        {
            "columns": {"smiles": {"dtype": "str"}, "inchikey": {"dtype": "str"}},
            "checks": [
                {"inchikey_matches": {"smiles_column": "smiles", "inchikey_column": "inchikey"}}
            ],
        }
    )
    table = pd.DataFrame({"smiles": ["CCO", "CCN"], "inchikey": [pc.inchikey("CCO")] * 2})
    assert pc.split(table, schema).quarantined["reasons"].to_dict() == {1: "inchikey_matches"}


def test_rules_are_read_from_yaml(tmp_path):
    path = tmp_path / "rules.yaml"
    path.write_text(EXAMPLE_RULES, encoding="utf-8")
    schema = pc.load_schema(path)
    assert set(schema.columns) == {"compound_id", "smiles", "pic50"}
    names = [c.name for c in schema.columns["smiles"].checks]
    assert names == [
        "valid_smiles",
        "allowed_elements",
        "single_fragment",
        "unique_molecules",
        "heavy_atoms",
        "stereo_defined",
    ]
    assert schema.columns["smiles"].checks[-1].raise_warning


def test_every_chemistry_check_can_be_named_in_a_rules_file():
    assert set(CHEM_CHECKS) == {n for n in pc.__all__ if n in CHEM_CHECKS}
    assert len(CHEM_CHECKS) == 11


def test_dtypes_and_column_settings():
    schema = pc.schema_from_dict(
        {
            "columns": {"n": {"dtype": "int", "coerce": True, "description": "a count"}},
            "strict": True,
        }
    )
    assert schema.strict
    assert schema.columns["n"].coerce
    assert schema.columns["n"].description == "a count"
    assert schema.validate(pd.DataFrame({"n": ["1", "2"]}))["n"].tolist() == [1, 2]


@pytest.mark.parametrize(
    ("rules", "message"),
    [
        ([], "must be a mapping"),
        ({"columns": {}}, "at least one column"),
        ({"columns": {"x": {}}, "colums": 1}, "unknown top-level"),
        ({"columns": {"x": {"dtyp": "str"}}}, "unknown settings"),
        ({"columns": {"x": {"dtype": "string"}}}, "dtype must be one of"),
        ({"columns": {"x": "str"}}, "must be a mapping"),
        ({"columns": {"x": {"checks": "valid_smiles"}}}, "must be a list"),
        ({"columns": {"x": {"checks": ["valid_smile"]}}}, "unknown check 'valid_smile'"),
        ({"columns": {"x": {"checks": [{"heavy_atoms": 5}]}}}, "must be a mapping"),
        ({"columns": {"x": {"checks": [{"heavy_atoms": {"min": 5}}]}}}, "unknown arguments"),
        ({"columns": {"x": {"checks": [{"in_range": {"min_value": 1}}]}}}, "in_range"),
        ({"columns": {"x": {"checks": [{"a": 1, "b": 2}]}}}, "a name with its arguments"),
        (
            {"columns": {"x": {"checks": [{"no_substructures": {"patterns": {"bad": "C(("}}}]}}},
            "SMARTS",
        ),
    ],
)
def test_bad_rules_say_where_and_why(rules, message):
    with pytest.raises(pc.RulesError, match=message):
        pc.schema_from_dict(rules)


def test_a_check_with_no_arguments_written_as_a_mapping():
    schema = pc.schema_from_dict({"columns": {"s": {"checks": [{"valid_smiles": None}]}}})
    assert schema.columns["s"].checks[0].name == "valid_smiles"


def test_invalid_yaml_is_reported(tmp_path):
    path = tmp_path / "rules.yaml"
    path.write_text("columns: [unclosed", encoding="utf-8")
    with pytest.raises(pc.RulesError, match="not valid YAML"):
        pc.load_schema(path)


def test_pandera_checks_are_allowed():
    schema = pc.schema_from_dict(
        {"columns": {"u": {"checks": [{"isin": {"allowed_values": ["nM", "uM"]}}]}}}
    )
    with pytest.raises(pa.errors.SchemaError):
        schema.validate(pd.DataFrame({"u": ["nM", "mM"]}))
