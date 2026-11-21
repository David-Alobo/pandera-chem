"""Every example in the docs is real: the tests run it and compare what the page shows.

Three kinds of marker, each just before a fenced block:

- ``<!-- example: examples/name.py -->``: the block shows the script (without its docstring); an
  optional ``text`` block after it shows what it prints.
- ``<!-- file: examples/name.ext -->``: the block shows the file exactly.
- ``<!-- command: pandera-chem ... -->``: the next ``text`` block shows what the command prints,
  run from a copy of the repository's ``examples`` folder.
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = [ROOT / "README.md", ROOT / "docs" / "usage.md"]
EXAMPLE = re.compile(
    r"<!-- example: (\S+) -->\s*```python\n(.*?)```"
    r"(?:\s*(?:It prints:)?\s*```text\n(.*?)```)?",
    re.S,
)
FILE = re.compile(r"<!-- file: (\S+) -->\s*```\w*\n(.*?)```", re.S)
COMMAND = re.compile(r"<!-- command: (.+?) -->.*?```text\n(.*?)```", re.S)


def found(pattern):
    """Every match of one marker in the docs, as pytest parameters named after the doc."""
    params = []
    for doc in DOCS:
        for m in pattern.finditer(doc.read_text(encoding="utf-8")):
            params.append(pytest.param(*m.groups(), id=f"{doc.name}:{m.group(1)}"))
    return params


def test_every_example_script_is_shown_somewhere():
    shown = {p.values[0] for p in found(EXAMPLE)}
    scripts = {f"examples/{p.name}" for p in (ROOT / "examples").glob("*.py")}
    assert scripts == shown


def test_the_docs_have_files_and_commands_to_check():
    assert found(FILE) and found(COMMAND)


@pytest.mark.parametrize(("script", "code", "output"), found(EXAMPLE))
def test_docs_show_the_script_and_what_it_prints(script, code, output, tmp_path):
    source = (ROOT / script).read_text(encoding="utf-8")
    assert code == source.split('"""', 2)[2].lstrip("\n"), "the code shown differs from the script"
    run = subprocess.run(
        [sys.executable, str(ROOT / script)], capture_output=True, text=True, cwd=tmp_path
    )
    assert run.returncode == 0, run.stderr
    if output is not None:
        assert output == run.stdout


@pytest.mark.parametrize(("path", "shown"), found(FILE))
def test_docs_show_the_file_exactly(path, shown):
    assert shown == (ROOT / path).read_text(encoding="utf-8")


@pytest.mark.parametrize(("command", "output"), found(COMMAND))
def test_docs_show_what_the_command_prints(command, output, tmp_path):
    shutil.copytree(ROOT / "examples", tmp_path / "examples")
    program, *args = command.split()
    assert program == "pandera-chem"
    run = subprocess.run(
        [sys.executable, "-m", "pandera_chem", *args], capture_output=True, text=True, cwd=tmp_path
    )
    assert run.returncode in (0, 1), run.stderr
    assert output == run.stdout
