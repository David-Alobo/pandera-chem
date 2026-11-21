import json

import pandas as pd
import pandera.pandas as pa

import pandera_chem as pc
from pandera_chem import report as report_module

SCHEMA = pa.DataFrameSchema(
    {"smiles": pa.Column(str, [pc.valid_smiles(), pc.stereo_defined(raise_warning=True)])}
)


def result_for(smiles):
    return pc.split(pd.DataFrame({"smiles": smiles}), SCHEMA)


def test_json_report_has_the_counts():
    data = json.loads(pc.report_json(result_for(["CCO", "c1ccccc", "CC(N)C(=O)O"]), "in.csv"))
    assert data["source"] == "in.csv"
    assert (data["rows_in"], data["rows_valid"], data["rows_quarantined"]) == (3, 2, 1)
    assert data["reconciles"] is True
    assert data["by_reason"] == {"smiles: valid_smiles": 1}
    assert len(data["warnings"]) == 1
    assert data["pandera_chem_version"] == pc.__version__


def test_html_report_shows_counts_reasons_and_rows():
    html = pc.report_html(result_for(["CCO", "c1ccccc"]), "My report", "in.csv")
    assert html.startswith("<!doctype html>")
    assert "<title>My report</title>" in html
    assert "<code>in.csv</code>" in html
    assert "smiles: valid_smiles" in html
    assert "c1ccccc" in html
    assert "quarantined (50.0%)" in html
    assert "<script" not in html  # self-contained and inert


def test_html_report_escapes_the_data():
    html = pc.report_html(result_for(["<script>alert(1)</script>"]), "<b>t</b>", "<x>.csv")
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert "<b>t</b>" not in html


def test_html_report_when_everything_passes():
    html = pc.report_html(result_for(["CCO"]))
    assert "No row failed any check." in html
    assert "No rows were quarantined." in html
    assert "No warnings." in html
    assert "Source: a data frame" in html


def test_html_report_shows_only_the_first_rows(monkeypatch):
    monkeypatch.setattr(report_module, "MAX_ROWS", 2)
    html = pc.report_html(result_for(["c1ccccc", "c1cccc", "C1CC", "CCO"]))
    assert "Showing the first 2 of 3 rows" in html
