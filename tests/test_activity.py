import importlib
import sys

import pandas as pd
import pandera.pandas as pa
import pytest

import pandera_chem as pc
from pandera_chem import _units

MU, MICRO = "μM", "µM"  # Greek mu and the micro sign: they look the same


def failing_rows(schema, df):
    """The index of every row a schema fails, as a sorted list."""
    try:
        schema.validate(df, lazy=True)
        return []
    except pa.errors.SchemaErrors as error:
        return sorted(set(error.failure_cases["index"].dropna().astype(int)))


@pytest.mark.parametrize(
    ("units", "ok"),
    [
        ("nM", True),
        ("uM", True),
        (MICRO, True),
        (MU, True),
        ("mM", True),
        ("pM", True),
        ("mol/L", True),
        ("umol/L", True),
        ("nm", False),  # a nanometre
        ("ug/mL", False),  # a mass per volume
        ("%", False),
        ("n/a", False),
        ("10 nM", False),
        ("", False),
    ],
)
def test_concentration_units(units, ok):
    schema = pa.DataFrameSchema({"u": pa.Column(str, pc.concentration_units())})
    assert failing_rows(schema, pd.DataFrame({"u": [units]})) == ([] if ok else [0])


def test_nm_per_unit_factors():
    assert _units.nm_per_unit("nM") == pytest.approx(1)
    assert _units.nm_per_unit(MICRO) == pytest.approx(1_000)
    assert _units.nm_per_unit("mM") == pytest.approx(1e6)
    assert _units.nm_per_unit("M") == pytest.approx(1e9)
    assert _units.nm_per_unit("pM") == pytest.approx(1e-3)
    assert _units.nm_per_unit(None) is None
    _units.clear_cache()
    assert _units.nm_per_unit("uM") == pytest.approx(1_000)


TABLE = pd.DataFrame(
    {
        "value": [50, 0.05, 4.2, 50, 1e-6, -3, 10, 250, 100],
        "units": ["nM", MICRO, "uM", "mM", "nM", "nM", "nm", "ug/mL", "nM"],
        "relation": ["=", "=", ">", "=", "=", "=", "=", "=", "="],
        "pchembl": [7.30, 7.30, 9.99, 7.30, None, None, 8.0, 6.0, None],
    }
)


def test_pchembl_matches_compares_only_exact_usable_rows():
    schema = pa.DataFrameSchema(
        checks=[pc.pchembl_matches("value", "units", "pchembl", relation_column="relation")]
    )
    # row 2 is ">" (not compared); rows 4-5 have no pChEMBL; rows 6-7 have no molar unit
    assert failing_rows(schema, TABLE) == [3]  # 50 mM, but a pChEMBL that fits 50 nM


def test_pchembl_matches_without_a_relation_column_compares_every_row():
    schema = pa.DataFrameSchema(checks=[pc.pchembl_matches("value", "units", "pchembl")])
    assert failing_rows(schema, TABLE) == [2, 3]  # now the ">" row is compared too


def test_pchembl_tolerance():
    df = pd.DataFrame({"value": [100.0], "units": ["nM"], "pchembl": [7.03]})
    strict = pa.DataFrameSchema(checks=[pc.pchembl_matches("value", "units", "pchembl")])
    loose = pa.DataFrameSchema(
        checks=[pc.pchembl_matches("value", "units", "pchembl", tolerance=0.05)]
    )
    assert failing_rows(strict, df) == [0]
    assert failing_rows(loose, df) == []


def test_activity_in_range():
    schema = pa.DataFrameSchema(checks=[pc.activity_in_range("value", "units")])
    assert failing_rows(schema, TABLE) == [4, 5]  # below 1 pM, and negative
    narrow = pa.DataFrameSchema(checks=[pc.activity_in_range("value", "units", max_nm=1e5)])
    assert failing_rows(narrow, TABLE) == [3, 4, 5]  # 50 mM is now too weak


def test_values_given_as_text_are_read_as_numbers():
    df = pd.DataFrame({"value": ["50", "n/a"], "units": ["nM", "nM"], "pchembl": ["7.30", ""]})
    schema = pa.DataFrameSchema(
        checks=[
            pc.pchembl_matches("value", "units", "pchembl"),
            pc.activity_in_range("value", "units"),
        ]
    )
    assert failing_rows(schema, df) == []


def test_each_problem_is_reported_once_in_split():
    schema = pa.DataFrameSchema(
        {
            "value": pa.Column(float),
            "units": pa.Column(str, pc.concentration_units()),
            "relation": pa.Column(str),
            "pchembl": pa.Column(float, nullable=True),
        },
        checks=[
            pc.pchembl_matches("value", "units", "pchembl", relation_column="relation"),
            pc.activity_in_range("value", "units"),
        ],
    )
    reasons = pc.split(TABLE, schema).quarantined["reasons"].to_dict()
    assert reasons == {
        3: "pchembl_matches",
        4: "activity_in_range",
        5: "activity_in_range",
        6: "units: concentration_units",
        7: "units: concentration_units",
    }


def test_activity_checks_from_a_rules_file():
    schema = pc.schema_from_dict(
        {
            "columns": {
                "value": {"dtype": "float"},
                "units": {"dtype": "str", "checks": ["concentration_units"]},
                "pchembl": {"dtype": "float", "nullable": True},
            },
            "checks": [
                {
                    "pchembl_matches": {
                        "value_column": "value",
                        "units_column": "units",
                        "pchembl_column": "pchembl",
                    }
                },
                {
                    "activity_in_range": {
                        "value_column": "value",
                        "units_column": "units",
                        "max_nm": 1e6,
                    }
                },
            ],
        }
    )
    df = pd.DataFrame({"value": [50.0, 50.0], "units": ["nM", "mM"], "pchembl": [7.30, None]})
    assert pc.split(df, schema).quarantined["reasons"].to_dict() == {1: "activity_in_range"}


def test_without_pint_the_checks_say_how_to_install_it(monkeypatch):
    _units.registry.cache_clear()
    monkeypatch.setitem(sys.modules, "pint", None)  # makes "import pint" fail
    try:
        with pytest.raises(ImportError, match=r"pandera-chem\[units\]"):
            pc.concentration_units()
        importlib.import_module("pandera_chem")  # the rest of the package still imports
        assert pc.valid_smiles().name == "valid_smiles"
    finally:
        monkeypatch.delitem(sys.modules, "pint")
        _units.registry.cache_clear()
