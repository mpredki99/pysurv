# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from typing import Any

from pysurveyor.backend.backend import (
    NA,
    NAType,
    Series,
    is_bool,
    is_number,
    isna,
    to_numeric,
)
from pysurveyor.validators.abc.conversion_validator import ConversionValidator


# ======================================================================================
# Numeric type validator
# ======================================================================================
class IsNumeric(ConversionValidator):
    """Validate and convert values to a numeric Series."""

    _output_dtype = int
    _fill_value = float("nan")

    # ==================================================================================
    def _convert(self, data: Series) -> Series:
        """Convert values to numbers and coerce invalid values to missing."""
        return to_numeric(data, errors="coerce")


# ======================================================================================
# Text type validator
# ======================================================================================
class IsText(ConversionValidator):
    """Validate and convert values to a text Series."""

    _output_dtype = "str"
    _fill_value = NA

    # ==================================================================================
    def _convert(self, data: Series) -> Series:
        """Convert values to strings."""
        return data.astype("str")


# ======================================================================================
# Boolean type validator
# ======================================================================================
class IsBoolean(ConversionValidator):
    """Validate and convert values to a boolean Series."""

    _output_dtype = "boolean"
    _fill_value = False

    # ==================================================================================
    def _convert(self, data: Series) -> Series:
        """Convert values to boolean."""

        # Define a helper function to convert individual values to boolean or NA -------
        def to_bool_or_na(value: Any) -> bool | NAType:
            if isna(value):
                return NA

            # --------------------------------------------------------------------------
            if is_bool(value):
                return bool(value)

            # --------------------------------------------------------------------------
            if is_number(value):
                return {0: False, 1: True}.get(value, NA)

            # --------------------------------------------------------------------------
            if isinstance(value, str):
                value = value.strip().lower()

                if value in {"true", "t"}:
                    return True
                if value in {"false", "f"}:
                    return False

                # ----------------------------------------------------------------------
                try:
                    number = float(value)
                except ValueError:
                    return NA

                # ----------------------------------------------------------------------
                if number == 1:
                    return True
                if number == 0:
                    return False

            # --------------------------------------------------------------------------
            return NA

        # ------------------------------------------------------------------------------
        return data.map(to_bool_or_na).astype("boolean")
