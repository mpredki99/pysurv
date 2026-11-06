# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from typing import Any

import pandas as pd

from pysurveyor.validators.abc.conversion_validator import ConversionValidator


# ======================================================================================
# Numeric type validator
# ======================================================================================
class IsNumeric(ConversionValidator):
    """Validate and convert values to a numeric Series."""

    _output_dtype = int
    _fill_value = float("nan")

    # ==================================================================================
    def _convert(self, data: pd.Series) -> pd.Series:
        """Convert values to numbers and coerce invalid values to missing."""
        return pd.to_numeric(data, errors="coerce")


# ======================================================================================
# Text type validator
# ======================================================================================
class IsText(ConversionValidator):
    """Validate and convert values to a text Series."""

    _output_dtype = "str"
    _fill_value = pd.NA

    # ==================================================================================
    def _convert(self, data: pd.Series) -> pd.Series:
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
    def _convert(self, data: pd.Series) -> pd.Series:
        """Convert values to boolean."""

        # Define a helper function to convert individual values to boolean or NA -------
        def to_bool_or_na(value: Any) -> bool | type(pd.NA):  # type: ignore
            if pd.isna(value):
                return pd.NA

            # --------------------------------------------------------------------------
            if isinstance(value, bool):
                return value

            # --------------------------------------------------------------------------
            if pd.api.types.is_number(value):
                if value == 0:
                    return False
                if value == 1:
                    return True
                return pd.NA

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
                    return pd.NA

                # ----------------------------------------------------------------------
                if number == 1:
                    return True
                if number == 0:
                    return False

            # --------------------------------------------------------------------------
            return pd.NA

        # ------------------------------------------------------------------------------
        return data.map(to_bool_or_na).astype("boolean")
