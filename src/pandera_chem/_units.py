"""Reading concentration units with pint, once per spelling.

pint is optional (``pip install "pandera-chem[units]"``): only the activity checks need it, and
they say how to install it if it is missing.
"""

from __future__ import annotations

import importlib
from functools import lru_cache
from typing import Any

INSTALL = 'the activity checks need pint: pip install "pandera-chem[units]"'


@lru_cache(maxsize=1)
def registry() -> Any:
    """The unit registry from pint, made once.

    Returns:
        The registry.

    Raises:
        ImportError: If pint is not installed, saying how to install it.
    """
    try:
        pint = importlib.import_module("pint")
    except ImportError as error:
        raise ImportError(INSTALL) from error
    return pint.UnitRegistry()


@lru_cache(maxsize=4096)
def nm_per_unit(units: str) -> float | None:
    """How many nanomolar one of a unit is, if the unit is a molar concentration.

    Args:
        units: A unit as written, e.g. "uM", "µM", "mol/L".

    Returns:
        The factor to nanomolar (1000.0 for micromolar), or None if pint cannot read the
        spelling or it is not an amount per volume (a mass per volume such as "ug/mL", a length
        such as "nm", or a percentage).
    """
    if not isinstance(units, str) or not units.strip():
        return None
    ureg = registry()
    try:
        unit = ureg.Unit(units)
        if unit.dimensionality != ureg.Unit("M").dimensionality:
            return None
        return float(ureg.Quantity(1, unit).to("nM").magnitude)
    except (ValueError, TypeError, AttributeError, KeyError):
        # pint's own errors (an unknown unit, a scaling factor in the unit) subclass these
        return None


def clear_cache() -> None:
    """Forget the units already read.

    Returns:
        None.
    """
    nm_per_unit.cache_clear()
