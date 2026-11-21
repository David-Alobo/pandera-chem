import tomllib
from pathlib import Path

import pandera_chem

ROOT = Path(__file__).resolve().parents[1]


def test_version_matches_pyproject():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    assert pandera_chem.__version__ == project["version"]
