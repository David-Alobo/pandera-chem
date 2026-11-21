"""Usage guide: the first step of a pipeline, which keeps good rows and saves the rest."""

import datetime
import logging
import sys
from pathlib import Path

import pandas as pd
import pandera.pandas as pa

import pandera_chem as pc

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s", stream=sys.stdout)
log = logging.getLogger("ingest")

SCHEMA = pa.DataFrameSchema(
    {
        "compound_id": pa.Column(str, unique=True),
        "smiles": pa.Column(str, [pc.valid_smiles(), pc.single_fragment(), pc.unique_molecules()]),
        "pic50": pa.Column(float, pa.Check.in_range(0, 14)),
    }
)
MAX_QUARANTINED = 0.25  # a larger share usually means something is wrong upstream


def ingest(path: Path, out: Path) -> pd.DataFrame:
    """Read a supplier's file, keep the rows that pass, and save the others with their reasons."""
    df = pd.read_csv(path).reset_index(drop=True)  # split() needs a unique index
    result = pc.split(df, SCHEMA)
    log.info("%s\n%s", path.name, result.report)
    if not result.report.reconciles:
        raise RuntimeError("rows went missing")
    stamp = datetime.date(2026, 11, 21).isoformat()  # in a real job: datetime.date.today()
    out.mkdir(exist_ok=True)
    result.valid.to_csv(out / f"valid-{stamp}.csv", index=False)
    result.quarantined.to_csv(out / f"quarantine-{stamp}.csv", index=False)
    share = result.report.rows_quarantined / result.report.rows_in
    if share > MAX_QUARANTINED:
        raise RuntimeError(f"{share:.0%} of rows quarantined; check the supplier's file")
    return result.valid


# a supplier's file, as it might arrive
pd.DataFrame(
    {
        "compound_id": ["S-1", "S-2", "S-3", "S-4", "S-5", "S-6", "S-7", "S-8"],
        "smiles": ["CCO", "CCN", "c1ccccc1O", "CCN.Cl", "CC(=O)O", "OCC", "CCCC", "CCOC"],
        "pic50": [5.1, 6.3, 4.8, 5.5, 7.0, 5.1, 4.2, 6.6],
    }
).to_csv("supplier.csv", index=False)

clean = ingest(Path("supplier.csv"), Path("out"))
print(len(clean), "rows go on to the next step")
print(sorted(p.name for p in Path("out").iterdir()))
