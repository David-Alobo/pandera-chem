"""Every example in the docs is a script in examples/, and the docs show its code and output."""

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = [ROOT / "README.md", ROOT / "docs" / "usage.md"]
# <!-- example: examples/name.py -->, then a ```python block, optionally followed by a ```text block
MARKED = re.compile(
    r"<!-- example: (\S+) -->\s*```python\n(.*?)```"
    r"(?:\s*(?:It prints:)?\s*```text\n(.*?)```)?",
    re.S,
)


def examples():
    """Each marked example: the doc, the script, the code shown and the output shown."""
    found = []
    for doc in DOCS:
        for m in MARKED.finditer(doc.read_text(encoding="utf-8")):
            found.append(
                pytest.param(m.group(1), m.group(2), m.group(3), id=f"{doc.name}:{m.group(1)}")
            )
    return found


def test_every_example_script_is_shown_somewhere():
    shown = {p.values[0] for p in examples()}
    scripts = {f"examples/{p.name}" for p in (ROOT / "examples").glob("*.py")}
    assert scripts == shown


@pytest.mark.parametrize(("script", "code", "output"), examples())
def test_docs_show_the_script_and_what_it_prints(script, code, output, tmp_path):
    source = (ROOT / script).read_text(encoding="utf-8")
    assert code == source.split('"""', 2)[2].lstrip("\n"), "the code shown differs from the script"
    run = subprocess.run(
        [sys.executable, str(ROOT / script)], capture_output=True, text=True, cwd=tmp_path
    )
    assert run.returncode == 0, run.stderr
    if output is not None:
        assert output == run.stdout
