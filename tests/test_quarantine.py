import warnings

import pandas as pd
import pandera.pandas as pa
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pandera_chem as pc

SCHEMA = pa.DataFrameSchema(
    {
        "smiles": pa.Column(
            str,
            [
                pc.valid_smiles(),
                pc.allowed_elements(),
                pc.single_fragment(),
                pc.neutral(),
                pc.unique_molecules(),
                pc.stereo_defined(raise_warning=True),
            ],
        ),
        "pic50": pa.Column(float, pa.Check.in_range(0, 14)),
    }
)
TABLE = pd.DataFrame(
    {
        "smiles": [
            "CCO",
            "OCC",
            "c1ccccc1",
            "c1ccccc",
            "C[Se]C",
            "CC(=O)[O-].[Na+]",
            "CC(N)C(=O)O",
        ],
        "pic50": [5.0, 6.0, 7.0, 4.0, 5.0, 16.0, 5.0],
    }
)


def test_split_keeps_good_rows_and_explains_bad_ones():
    result = pc.split(TABLE, SCHEMA)
    assert list(result.valid.index) == [0, 2, 6]
    assert result.quarantined["reasons"].to_dict() == {
        1: "smiles: unique_molecules",
        3: "smiles: valid_smiles",
        4: "smiles: allowed_elements",
        5: "smiles: allowed_elements; smiles: single_fragment; pic50: in_range(0, 14)",
    }


def test_report_counts_reconcile():
    report = pc.split(TABLE, SCHEMA).report
    assert (report.rows_in, report.rows_valid, report.rows_quarantined) == (7, 3, 4)
    assert report.reconciles
    assert report.by_reason["smiles: allowed_elements"] == 2
    assert "7 rows in: 3 valid, 4 quarantined" in str(report)


def test_warnings_are_reported_not_raised_or_quarantined():
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning that escapes split() fails the test
        result = pc.split(TABLE, SCHEMA)
    assert 6 in result.valid.index
    assert len(result.report.warnings) == 1
    assert "stereo_defined" in result.report.warnings[0]


def test_quarantined_rows_are_unchanged():
    result = pc.split(TABLE, SCHEMA)
    pd.testing.assert_frame_equal(
        result.quarantined.drop(columns="reasons"), TABLE.loc[[1, 3, 4, 5]]
    )


def test_a_clean_table_is_all_valid():
    clean = TABLE.loc[[0, 2]]
    result = pc.split(clean, SCHEMA)
    assert len(result.valid) == 2
    assert result.quarantined.empty
    assert "reasons" in result.quarantined.columns
    assert result.report.reconciles


def test_a_table_wide_failure_raises():
    with pytest.raises(pa.errors.SchemaErrors):
        pc.split(TABLE.drop(columns="pic50"), SCHEMA)


def test_a_repeated_index_is_refused():
    with pytest.raises(ValueError, match="repeated labels"):
        pc.split(TABLE.set_axis([0, 0, 1, 2, 3, 4, 5]), SCHEMA)


def test_a_non_default_index_is_kept():
    table = TABLE.set_axis([f"ML-{i:03d}" for i in range(len(TABLE))])
    result = pc.split(table, SCHEMA)
    assert list(result.quarantined.index) == ["ML-001", "ML-003", "ML-004", "ML-005"]


POOL = ["CCO", "OCC", "c1ccccc1", "c1ccccc", "C[Se]C", "CC[NH3+]", "CCN", "CC(=O)O", "C1CC1"]


@settings(max_examples=50, deadline=None)
@given(st.lists(st.sampled_from(POOL), min_size=1, max_size=30))
def test_split_is_a_partition(picks):
    table = pd.DataFrame({"smiles": picks, "pic50": 5.0})
    result = pc.split(table, SCHEMA)
    kept, set_aside = set(result.valid.index), set(result.quarantined.index)
    assert kept.isdisjoint(set_aside)
    assert kept | set_aside == set(table.index)
    assert result.report.reconciles
