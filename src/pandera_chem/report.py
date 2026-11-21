"""Report cards: what happened to a table, as a web page for people and JSON for machines.

The HTML report is one self-contained file (no scripts, no downloads), so it can be attached to a
CI run, emailed or opened offline. Everything taken from the data is escaped.
"""

from __future__ import annotations

import json
from html import escape
from typing import Any

from . import __version__
from .quarantine import REASONS, Split

MAX_ROWS = 200  # quarantined rows shown in the HTML report; the CSV has them all

CSS = """
:root { --bg: #ffffff; --fg: #1d2030; --muted: #5c6178; --line: #e3e5ee; --card: #f6f7fb;
        --ok: #0f7b4f; --bad: #b42318; --warn: #a15c07; --bar: #7c5cff; }
@media (prefers-color-scheme: dark) {
  :root { --bg: #0d0f16; --fg: #eef0f5; --muted: #9aa0b8; --line: #262a3b; --card: #151826;
          --ok: #34d399; --bad: #f87171; --warn: #fbbf24; --bar: #a996ff; }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--fg);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1000px; margin: 0 auto; padding: 32px 16px 48px; }
h1 { font-size: 26px; margin: 0 0 4px; } h2 { font-size: 18px; margin: 32px 0 12px; }
.meta { color: var(--muted); margin: 0 0 24px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; }
.tile { background: var(--card); border: 1px solid var(--line); border-radius: 12px;
        padding: 14px; }
.tile b { display: block; font-size: 28px; font-variant-numeric: tabular-nums; }
.tile span { color: var(--muted); }
.ok { color: var(--ok); } .bad { color: var(--bad); } .warn { color: var(--warn); }
table { width: 100%; border-collapse: collapse; font-size: 14px; }
th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid var(--line);
         vertical-align: top; }
th { color: var(--muted); font-weight: 600; }
td.n { text-align: right; font-variant-numeric: tabular-nums; width: 5em; }
.reasons td { vertical-align: middle; }
.bar { height: 8px; border-radius: 4px; background: var(--bar); min-width: 2px; }
.scroll { overflow-x: auto; border: 1px solid var(--line); border-radius: 12px; }
.scroll table td, .scroll table th { white-space: nowrap; }
code { font: 13px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
ul.warnings { padding-left: 20px; } ul.warnings li { margin-bottom: 6px; }
footer { margin-top: 40px; color: var(--muted); font-size: 13px; }
"""


def report_dict(result: Split, source: str | None = None) -> dict[str, Any]:
    """The report as plain data, for JSON and CI systems.

    Args:
        result: What :func:`pandera_chem.split` returned.
        source: Where the table came from (e.g. a file name), if known.

    Returns:
        The counts, the failures by reason, the warnings and whether the counts reconcile.
    """
    r = result.report
    return {
        "source": source,
        "rows_in": r.rows_in,
        "rows_valid": r.rows_valid,
        "rows_quarantined": r.rows_quarantined,
        "reconciles": r.reconciles,
        "by_reason": r.by_reason,
        "warnings": r.warnings,
        "pandera_chem_version": __version__,
    }


def report_json(result: Split, source: str | None = None) -> str:
    """The report as JSON text.

    Args:
        result: What :func:`pandera_chem.split` returned.
        source: Where the table came from, if known.

    Returns:
        Indented JSON.
    """
    return json.dumps(report_dict(result, source), indent=2, ensure_ascii=False) + "\n"


def _tile(value: object, label: str, css: str = "") -> str:
    """One summary tile.

    Args:
        value: The number or word to show large.
        label: What it is.
        css: An extra class for its colour.

    Returns:
        The tile's HTML.
    """
    shown = f'<b class="{css}">{escape(str(value))}</b>'
    return f'<div class="tile">{shown}<span>{escape(label)}</span></div>'


def _reasons_table(by_reason: dict[str, int], rows_in: int) -> str:
    """The table of reasons, each with its count and a bar.

    Args:
        by_reason: How many rows failed each reason.
        rows_in: Rows in the table, for the bars' scale.

    Returns:
        The table's HTML, or a sentence if nothing failed.
    """
    if not by_reason:
        return "<p>No row failed any check.</p>"
    rows = []
    for reason, n in by_reason.items():
        width = 100 * n / max(rows_in, 1)
        rows.append(
            f'<tr><td><code>{escape(reason)}</code></td><td class="n">{n}</td>'
            f'<td style="width:40%"><div class="bar" style="width:{width:.1f}%"></div></td></tr>'
        )
    head = "<tr><th>Reason</th><th class='n'>Rows</th><th>Share of all rows</th></tr>"
    return f"<table class='reasons'>{head}{''.join(rows)}</table>"


def _rows_table(result: Split) -> str:
    """The quarantined rows (the first :data:`MAX_ROWS`), reasons first.

    Args:
        result: What :func:`pandera_chem.split` returned.

    Returns:
        The table's HTML, or a sentence if no row was quarantined.
    """
    q = result.quarantined
    if q.empty:
        return "<p>No rows were quarantined.</p>"
    columns = [REASONS] + [c for c in q.columns if c != REASONS]
    head = "<tr><th>Row</th>" + "".join(f"<th>{escape(str(c))}</th>" for c in columns) + "</tr>"
    body = []
    for index, row in q.head(MAX_ROWS).iterrows():
        cells = "".join(f"<td>{escape(str(row[c]))}</td>" for c in columns)
        body.append(f"<tr><td>{escape(str(index))}</td>{cells}</tr>")
    more = ""
    if len(q) > MAX_ROWS:
        more = (
            f"<p>Showing the first {MAX_ROWS} of {len(q)} rows; "
            "the quarantine file has them all.</p>"
        )
    return f'<div class="scroll"><table>{head}{"".join(body)}</table></div>{more}'


def report_html(
    result: Split, title: str = "Data quality report", source: str | None = None
) -> str:
    """The report as one self-contained web page.

    Args:
        result: What :func:`pandera_chem.split` returned.
        title: The page's title.
        source: Where the table came from (e.g. a file name), shown under the title.

    Returns:
        The HTML document.
    """
    r = result.report
    share = r.rows_quarantined / r.rows_in if r.rows_in else 0.0
    tiles = "".join(
        [
            _tile(r.rows_in, "rows in"),
            _tile(r.rows_valid, "valid", "ok"),
            _tile(
                r.rows_quarantined,
                f"quarantined ({share:.1%})",
                "bad" if r.rows_quarantined else "",
            ),
            _tile(len(r.warnings), "warnings", "warn" if r.warnings else ""),
            _tile(
                "yes" if r.reconciles else "NO",
                "counts reconcile",
                "ok" if r.reconciles else "bad",
            ),
        ]
    )
    warnings = (
        "<ul class='warnings'>" + "".join(f"<li>{escape(w)}</li>" for w in r.warnings) + "</ul>"
        if r.warnings
        else "<p>No warnings.</p>"
    )
    meta = f"Source: <code>{escape(source)}</code>" if source else "Source: a data frame"
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<style>{CSS}</style>
</head>
<body>
<main>
<h1>{escape(title)}</h1>
<p class="meta">{meta}</p>
<div class="tiles">{tiles}</div>
<h2>Why rows were quarantined</h2>
{_reasons_table(r.by_reason, r.rows_in)}
<p class="meta">A row that fails two checks is counted under both.</p>
<h2>Warnings</h2>
{warnings}
<h2>Quarantined rows</h2>
{_rows_table(result)}
<footer>Made by pandera-chem {escape(__version__)}.</footer>
</main>
</body>
</html>
"""
