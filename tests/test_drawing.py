import pandas as pd
import pandera.pandas as pa
import pytest

import pandera_chem as pc
from pandera_chem.drawing import explain, first_rows, picture
from pandera_chem.report import smiles_column


@pytest.mark.parametrize(
    ("smiles", "failed", "stats", "atoms", "caption"),
    [
        ("C[Se]c1ccccc1", ["allowed_elements"], {}, {1}, "element not allowed: Se"),
        (
            "CCO",
            ["allowed_elements"],
            {"allowed_elements": {"elements": ["C", "H"]}},
            {2},
            "element not allowed: O",
        ),
        ("CCN(CC)CC.Cl", ["single_fragment"], {}, {7}, "extra piece: Cl"),
        ("CC[NH3+]", ["neutral"], {}, {2}, "net charge +1"),
        ("CC(N)C(=O)O", ["stereo_defined"], {}, {1}, "no configuration: 1 stereocentre"),
        (
            "O=Cc1ccc(N)cc1",
            ["no_substructures"],
            {"no_substructures": {"patterns": {"aldehyde": "[CX3H1](=O)[#6]"}}},
            {0, 1, 2},
            "contains: aldehyde",
        ),
        ("C" * 60, ["heavy_atoms"], {}, set(), "60 heavy atoms"),
        ("CCO.[Na+]", ["mol_weight"], {}, set(), "46 g/mol"),
    ],
)
def test_explain_blames_the_right_atoms(smiles, failed, stats, atoms, caption):
    mol, why = explain(smiles, failed, stats)
    assert mol is not None
    assert why.atoms == atoms
    assert why.captions == [caption]


def test_a_double_bond_without_configuration_is_a_bond():
    _, why = explain("CC=CC", ["stereo_defined"])
    assert why.bonds == {1}
    assert why.captions == ["no configuration: 1 double bond"]


def test_a_pattern_that_does_not_match_still_has_a_caption():
    _, why = explain("CCO", ["no_substructures"])
    assert why.captions == ["contains: a forbidden group"]


def test_unreadable_smiles_get_no_molecule_and_rdkit_s_reason():
    mol, why = explain("C(C)(C)(C)(C)C", ["valid_smiles"])
    assert mol is None
    assert why.captions[0].startswith("cannot be read: Explicit valence")


def test_duplicates_say_which_row_they_repeat():
    _, why = explain("OCC", ["unique_molecules"], repeats="C1")
    assert why.captions == ["same molecule as row C1"]
    _, why = explain("OCC", ["unique_molecules"])
    assert why.captions == ["same molecule as an earlier row"]


def test_picture_is_inline_svg_and_does_not_change_the_molecule():
    mol, why = explain("C[Se]c1ccccc1", ["allowed_elements"])
    before = mol.ToBinary()
    svg = picture(mol, why)
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert "<?xml" not in svg and "<script" not in svg
    assert mol.ToBinary() == before


def test_first_rows():
    assert first_rows({0: "CCO", 1: "OCC", 2: "CCN", 3: None}) == {
        pc.inchikey("CCO"): 0,
        pc.inchikey("CCN"): 2,
    }


SCHEMA = pa.DataFrameSchema(
    {
        "id": pa.Column(str),
        "structure": pa.Column(
            str, [pc.valid_smiles(), pc.allowed_elements(), pc.unique_molecules()]
        ),
        "pic50": pa.Column(float, pa.Check.in_range(0, 14)),
    }
)
TABLE = pd.DataFrame(
    {
        "id": ["A", "B", "C", "D", "E"],
        "structure": ["CCO", "OCC", "c1ccccc", "C[Se]C", "CCN"],
        "pic50": [5.0, 5.0, 5.0, 5.0, 20.0],
    }
)


def test_the_report_card_draws_the_quarantined_molecules():
    result = pc.split(TABLE, SCHEMA)
    assert result.schema is SCHEMA
    assert smiles_column(result) == "structure"
    html = pc.report_html(result)
    assert html.count("class='card'") == 4  # B, C, D and E
    assert html.count("<svg") == 3  # C cannot be drawn
    assert "cannot be drawn" in html
    assert "same molecule as row 0" in html
    assert "element not allowed: Se" in html
    assert "cannot be read: the SMILES could not be read" in html


def test_pictures_can_be_turned_off():
    html = pc.report_html(pc.split(TABLE, SCHEMA), pictures=False)
    assert "<svg" not in html and "class='pic'" not in html
    assert html.count("class='card'") == 4


def test_no_pictures_without_a_smiles_column():
    schema = pa.DataFrameSchema({"x": pa.Column(int, pa.Check.gt(0))})
    html = pc.report_html(pc.split(pd.DataFrame({"x": [1, -1]}), schema))
    assert "<svg" not in html


def test_a_column_called_smiles_is_drawn_even_without_molecule_checks():
    schema = pa.DataFrameSchema({"smiles": pa.Column(str), "n": pa.Column(int, pa.Check.gt(0))})
    result = pc.split(pd.DataFrame({"smiles": ["CCO", None], "n": [-1, -1]}), schema)
    assert smiles_column(result) == "smiles"
    html = pc.report_html(result)
    assert html.count("<svg") == 1  # the empty SMILES gets an empty cell
