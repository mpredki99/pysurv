# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

import operator

from pysurveyor.validators.abc.comparison_validator import ComparisonValidator


# ======================================================================================
# Comparison value validators
# ======================================================================================
class Equal(ComparisonValidator):
    """Validate values equal to a threshold."""

    _operator = operator.eq
    _operator_error_str = "!="


# ======================================================================================
class GreaterOrEqual(ComparisonValidator):
    """Validate values greater than or equal to a threshold."""

    _operator = operator.ge
    _operator_error_str = "<"


# ======================================================================================
class Greater(ComparisonValidator):
    """Validate values greater than a threshold."""

    _operator = operator.gt
    _operator_error_str = "<="


# ======================================================================================
class LessOrEqual(ComparisonValidator):
    """Validate values less than or equal to a threshold."""

    _operator = operator.le
    _operator_error_str = ">"


# ======================================================================================
class Less(ComparisonValidator):
    """Validate values less than a threshold."""

    _operator = operator.lt
    _operator_error_str = ">="
