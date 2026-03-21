# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from copy import deepcopy
from functools import lru_cache
from typing import Any
from warnings import warn

from pysurv.typing.typing import AngleUnit, DistanceUnit, PySurvUnit
from pysurv.validators.between import Between
from pysurv.validators.constants import COMMENT
from pysurv.validators.equal import Equal
from pysurv.validators.greater import Greater
from pysurv.validators.greater_equal import GreaterOrEqual
from pysurv.validators.is_numeric import IsNumeric
from pysurv.validators.is_text import IsText
from pysurv.validators.less import Less
from pysurv.validators.less_equal import LessOrEqual
from pysurv.validators.pysurv_validator import PySurvValidator

# --------------------------------------------------------------------------------------
# Create namespace for eval function
# --------------------------------------------------------------------------------------
NAMESPACE = {
    "Between": Between,
    "Equal": Equal,
    "Greater": Greater,
    "GreaterOrEqual": GreaterOrEqual,
    "IsNumeric": IsNumeric,
    "IsText": IsText,
    "Less": Less,
    "LessOrEqual": LessOrEqual,
    # Other keywords
    "COMMENT": COMMENT,
}


# --------------------------------------------------------------------------------------
def compile_validator(value: Any) -> Any:
    if isinstance(value, PySurvValidator):
        return value

    if isinstance(value, str) and value.strip():
        try:
            # Try to build PySurvValidator
            return deepcopy(_compile_validator(value))
        except (NameError, TypeError, SyntaxError):
            warn(f"Could not parse PySurvValidator from value: {value}")

    return value


# --------------------------------------------------------------------------------------
def compile_unit(value: Any) -> Any:
    if isinstance(value, PySurvUnit):
        return value

    try:
        # Try to build PySurvUnit
        return _compile_unit(value)
    except (ValueError, KeyError):
        warn(f"Could not parse PySurvUnit from value: {value}")

    return value


# --------------------------------------------------------------------------------------
# Caching repetitive values
# --------------------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _compile_validator(expr: str) -> PySurvValidator:
    """Compile and cache a PySurvValidator object from a validator expression string."""
    if not expr:
        raise ValueError("validator expression must be non-empty")
    return eval(expr, {}, NAMESPACE)


def _compile_unit(value: PySurvUnit | str) -> PySurvUnit:
    """Compile and cache a PySurvUnit object from a expression string."""
    try:
        return DistanceUnit(value)
    except ValueError:
        return AngleUnit(value)
