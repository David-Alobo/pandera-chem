"""Usage guide: rules from a YAML file, and a report card, from Python."""

from pathlib import Path

import pandas as pd

import pandera_chem as pc

HERE = Path(__file__).parent
schema = pc.load_schema(HERE / "rules.yaml")
table = pd.read_csv(HERE / "compounds.csv")

result = pc.split(table, schema)
Path("report.html").write_text(pc.report_html(result, source="compounds.csv"), encoding="utf-8")
Path("report.json").write_text(pc.report_json(result, source="compounds.csv"), encoding="utf-8")

print(result.report.rows_valid, "valid rows;", result.report.rows_quarantined, "quarantined")
print(pc.report_dict(result)["by_reason"])
