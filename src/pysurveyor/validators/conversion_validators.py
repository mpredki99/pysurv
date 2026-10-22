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

    # ==================================================================================
    # Conversion configuration
    # ==================================================================================
    @property
    def _output_dtype(self) -> type:
        """Return the dtype used for empty converted Series."""
        return int

    @property
    def _fill_value(self) -> float:
        """Return the missing value used for ignored entries."""
        return float("nan")

    # ==================================================================================
    def _convert(self, data: Series) -> Series:
        """Convert values to numbers and coerce invalid values to missing."""
        return to_numeric(data, errors="coerce")


# ======================================================================================
# Text type validator
# ======================================================================================
class IsText(ConversionValidator):
    """Validate and convert values to a text Series."""

    # ==================================================================================
    # Conversion configuration
    # ==================================================================================
    @property
    def _output_dtype(self) -> str:
        """Return the dtype used for empty converted Series."""
        return "str"

    @property
    def _fill_value(self) -> object:
        """Return the missing value used for ignored entries."""
        return NA

    # ==================================================================================
    def _convert(self, data: Series) -> Series:
        """Convert values to strings."""
        return data.astype("str")
