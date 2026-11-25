"""Checks for activity data: the measurements that travel with the molecules.

A bioactivity table (IC50, Ki, EC50...) breaks in its own ways: a unit that is not a
concentration at all, a value off by a factor of 1,000 because of its unit, a stored pChEMBL that
disagrees with the value it came from. These checks need pint (``pip install
"pandera-chem[units]"``)::

    schema = pa.DataFrameSchema(
        {"value": pa.Column(float), "units": pa.Column(str, pc.concentration_units()),
         "pchembl": pa.Column(float, nullable=True)},
        checks=[pc.pchembl_matches("value", "units", "pchembl"),
                pc.activity_in_range("value", "units")],
    )

As with the molecule checks, each problem is reported once: a row whose unit cannot be read is
reported by ``concentration_units`` and passed by the table-wide checks.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pandera.pandas as pa

from ._units import nm_per_unit, registry


def concentration_units(**options: Any) -> pa.Check:
    """The unit is a molar concentration (an amount per volume), however it is spelled.

    Passes nM, uM, µM, μM, mM, mol/L, umol/L...; fails units that are not concentrations at all:
    a mass per volume (ug/mL, which needs a molecular weight to convert), a length (nm, a common
    slip for nM), a percentage, or text pint cannot read. To allow only certain spellings as
    well, add Pandera's ``isin`` check.

    Args:
        **options: Passed on to ``pandera.Check``, e.g. ``raise_warning=True`` to warn
            instead of failing, or ``n_failure_cases=10``.

    Returns:
        The check.
    """
    registry()  # fail now, with how to install pint, rather than on the first table
    return pa.Check(
        lambda units: nm_per_unit(units) is not None,
        element_wise=True,
        name="concentration_units",
        description="not a molar concentration unit",
        **options,
    )


def _in_nm(df: pd.DataFrame, value_column: str, units_column: str) -> pd.Series:
    """Each row's value in nanomolar, or NaN where the value or the unit cannot be used.

    Args:
        df: The table.
        value_column: The column of values.
        units_column: The column of units.

    Returns:
        The values in nanomolar.
    """
    values = pd.to_numeric(df[value_column], errors="coerce")
    factors = df[units_column].map(lambda u: nm_per_unit(u) if isinstance(u, str) else None)
    return values * pd.to_numeric(factors, errors="coerce")


def pchembl_matches(
    value_column: str,
    units_column: str,
    pchembl_column: str,
    relation_column: str | None = None,
    tolerance: float = 0.01,
    **options: Any,
) -> pa.Check:
    """A table-wide check: the stored pChEMBL equals 9 - log10(value in nM).

    pChEMBL is the negative log of the molar activity (pIC50 for an IC50), so 1 nM is 9 and
    1 µM is 6. Rows are compared only when every part is there and usable: a pChEMBL, a positive
    value, a molar unit and, if ``relation_column`` is given, the relation "=" (a value given as
    "> 10 µM" has no exact pChEMBL). The other rows pass; other checks report their problems.

    Args:
        value_column: The column of activity values.
        units_column: The column of their units.
        pchembl_column: The column of stored pChEMBL values.
        relation_column: The column of relations ("=", "<", ">"...), if the table has one.
        tolerance: How far the stored pChEMBL may be from the computed one.
        **options: Passed on to ``pandera.Check``.

    Returns:
        The check, to put in ``DataFrameSchema(checks=[...])``.
    """
    registry()

    def ok(df: pd.DataFrame) -> pd.Series:
        """Which rows' pChEMBL agree with their value and unit.

        Args:
            df: The table.

        Returns:
            True where they agree, or where they cannot be compared.
        """
        nm = _in_nm(df, value_column, units_column)
        stored = pd.to_numeric(df[pchembl_column], errors="coerce")
        comparable = stored.notna() & nm.notna() & (nm > 0)
        if relation_column is not None:
            comparable &= df[relation_column].astype(str).str.strip() == "="
        expected = 9 - np.log10(nm.where(comparable, 1.0))
        return ~comparable | ((stored - expected).abs() <= tolerance)

    return pa.Check(
        ok,
        name="pchembl_matches",
        description=f"{pchembl_column} is not 9 - log10({value_column} in nM)",
        statistics={"tolerance": tolerance},
        **options,
    )


def activity_in_range(
    value_column: str,
    units_column: str,
    min_nm: float = 0.001,
    max_nm: float = 1e9,
    **options: Any,
) -> pa.Check:
    """A table-wide check: each value, converted to nanomolar, is a plausible activity.

    The defaults allow 1 pM (0.001 nM) to 1 M (1e9 nM). A value outside that range is usually a
    unit slip: 50 entered with "mM" instead of "nM" is a million times too weak. A value of zero
    or below fails too. Rows whose value or unit cannot be used pass; other checks report them.

    Args:
        value_column: The column of activity values.
        units_column: The column of their units.
        min_nm: The smallest plausible value, in nanomolar.
        max_nm: The largest plausible value, in nanomolar.
        **options: Passed on to ``pandera.Check``.

    Returns:
        The check, to put in ``DataFrameSchema(checks=[...])``.
    """
    registry()

    def ok(df: pd.DataFrame) -> pd.Series:
        """Which rows' values are within the plausible range once in nanomolar.

        Args:
            df: The table.

        Returns:
            True where the value is plausible, or where it cannot be converted.
        """
        nm = _in_nm(df, value_column, units_column)
        return nm.isna() | ((nm >= min_nm) & (nm <= max_nm))

    return pa.Check(
        ok,
        name="activity_in_range",
        description=f"{value_column} outside {min_nm} to {max_nm} nM once converted",
        statistics={"min_nm": min_nm, "max_nm": max_nm},
        **options,
    )
