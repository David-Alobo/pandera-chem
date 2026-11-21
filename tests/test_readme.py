import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text(encoding="utf-8")
EXAMPLE = ROOT / "examples" / "quickstart.py"


def block(kind: str) -> str:
    """The first fenced block of one kind in the README."""
    return re.search(rf"```{kind}\n(.*?)```", README, re.S).group(1)


def test_readme_shows_the_example_code():
    code = EXAMPLE.read_text(encoding="utf-8")
    body = code.split('"""', 2)[2].lstrip("\n")  # without the module docstring
    assert block("python") == body


def test_readme_shows_what_the_example_prints():
    run = subprocess.run(
        [sys.executable, str(EXAMPLE)], capture_output=True, text=True, check=True, cwd=ROOT
    )
    assert block("text") == run.stdout


def test_version_matches_pyproject():
    import tomllib

    import pandera_chem

    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    assert pandera_chem.__version__ == project["version"]
