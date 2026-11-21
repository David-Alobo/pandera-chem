import json

import pandas as pd
import pytest
from rdkit import Chem

from pandera_chem import cli
from pandera_chem.config import EXAMPLE_RULES

CSV = (
    "compound_id,smiles,pic50\n"
    "C1,CC(=O)Oc1ccccc1C(=O)O,5.2\n"
    "C2,c1ccc2ccccc2c1,4.1\n"
    "C3,c1ccccc,6.0\n"
    "C4,CC(C)Cc1ccc(cc1)C(C)C(=O)O,n/a\n"
)


@pytest.fixture
def files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "rules.yaml").write_text(EXAMPLE_RULES, encoding="utf-8")
    (tmp_path / "compounds.csv").write_text(CSV, encoding="utf-8")
    return tmp_path


def test_check_writes_results_and_fails_on_any_quarantine_by_default(files, capsys):
    code = cli.main(["check", "compounds.csv", "--rules", "rules.yaml", "--out", "out"])
    printed = capsys.readouterr().out
    assert code == 1
    assert "4 rows in: 3 valid, 1 quarantined" in printed
    assert "FAIL: 25.0% of rows quarantined (allowed: 0.0%)" in printed
    out = files / "out"
    assert sorted(p.name for p in out.iterdir()) == [
        "compounds.quarantine.csv",
        "compounds.report.html",
        "compounds.report.json",
        "compounds.valid.csv",
    ]
    quarantine = pd.read_csv(out / "compounds.quarantine.csv")
    assert quarantine["reasons"].tolist() == ["smiles: valid_smiles"]
    assert json.loads((out / "compounds.report.json").read_text(encoding="utf-8"))["reconciles"]


def test_check_passes_within_the_allowed_share(files, capsys):
    code = cli.main(
        [
            "check",
            "compounds.csv",
            "--rules",
            "rules.yaml",
            "--out",
            "out",
            "--max-quarantined",
            "0.25",
        ]
    )
    assert code == 0
    assert "PASS: 25.0% of rows quarantined (allowed: 25.0%)" in capsys.readouterr().out


def test_text_that_is_not_a_number_becomes_empty_not_an_error(files):
    cli.main(["check", "compounds.csv", "--rules", "rules.yaml", "--out", "out"])
    valid = pd.read_csv(files / "out" / "compounds.valid.csv")
    assert valid.loc[valid["compound_id"] == "C4", "pic50"].isna().all()  # "n/a", nullable: true


def test_unusable_rules_or_data_give_exit_code_2(files, capsys):
    (files / "bad.yaml").write_text(
        "columns:\n  smiles:\n    checks: [valid_smile]\n", encoding="utf-8"
    )
    assert cli.main(["check", "compounds.csv", "--rules", "bad.yaml"]) == 2
    assert "unknown check 'valid_smile'" in capsys.readouterr().err
    assert cli.main(["check", "missing.csv", "--rules", "rules.yaml"]) == 2
    (files / "data.xlsx").write_text("x", encoding="utf-8")
    assert cli.main(["check", "data.xlsx", "--rules", "rules.yaml"]) == 2
    assert "cannot read .xlsx" in capsys.readouterr().err


def test_a_missing_column_gives_exit_code_2(files, capsys):
    (files / "nosmiles.csv").write_text("compound_id,pic50\nC1,5\n", encoding="utf-8")
    assert cli.main(["check", "nosmiles.csv", "--rules", "rules.yaml"]) == 2
    assert "cannot be checked as a whole" in capsys.readouterr().err


def test_tab_separated_files(files):
    (files / "compounds.tsv").write_text(CSV.replace(",", "\t"), encoding="utf-8")
    assert cli.main(["check", "compounds.tsv", "--rules", "rules.yaml", "--out", "out"]) == 1


def test_an_sdf_becomes_one_row_per_record(files):
    writer = Chem.SDWriter(str(files / "lib.sdf"))
    for i, smiles in enumerate(["CC(=O)Oc1ccccc1C(=O)O", "c1ccc2ccccc2c1"]):
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("compound_id", f"S{i}")
        mol.SetProp("pic50", "5.0")
        writer.write(mol)
    writer.close()
    with open(files / "lib.sdf", "a", encoding="utf-8") as f:  # a record RDKit cannot read
        f.write(
            "broken\n  RDKit\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0 0.0 0.0 Xx  0\n"
            "M  END\n> <compound_id>\nS9\n\n$$$$\n"
        )
    table = cli.read_table(files / "lib.sdf")
    assert table["compound_id"].tolist() == ["S0", "S1", "S9"]
    assert table["smiles"].tolist()[2] == ""
    assert cli.main(["check", "lib.sdf", "--rules", "rules.yaml", "--out", "out"]) == 1


def test_parquet_needs_pyarrow_or_reads(files):
    pd.DataFrame({"compound_id": ["C1"], "smiles": ["CCCCCCO"], "pic50": [5.0]}).pipe(
        lambda df: df.to_csv(files / "one.csv", index=False)
    )
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        (files / "one.parquet").write_bytes(b"PAR1")
        with pytest.raises(cli.InputError, match="pyarrow"):
            cli.read_table(files / "one.parquet")
    else:  # pragma: no cover - only when pyarrow is installed
        pd.read_csv(files / "one.csv").to_parquet(files / "one.parquet")
        assert len(cli.read_table(files / "one.parquet")) == 1


def test_init_writes_a_starter_file_and_will_not_overwrite(files, capsys):
    target = files / "new.yaml"
    assert cli.main(["init", str(target)]) == 0
    assert target.read_text(encoding="utf-8") == EXAMPLE_RULES
    assert cli.main(["init", str(target)]) == 2
    assert "already exists" in capsys.readouterr().err
    assert cli.main(["init", str(target), "--force"]) == 0


def test_version(capsys):
    with pytest.raises(SystemExit):
        cli.main(["--version"])
    assert "pandera-chem" in capsys.readouterr().out
