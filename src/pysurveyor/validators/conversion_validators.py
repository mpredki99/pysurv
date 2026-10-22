# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.backend.backend import NA, Series, to_numeric
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
