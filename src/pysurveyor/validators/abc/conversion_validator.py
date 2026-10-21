# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from abc import abstractmethod
from typing import Any

from pysurveyor.backend.backend import Series
from pysurveyor.validators.abc.validator import Validator


# ======================================================================================
# Conversion validator base class
# ======================================================================================
class ConversionValidator(Validator):
    """Base class for validators that validate values and converting their type.

    Concrete validators provide the output dtype, the value used for ignored
    entries, and the conversion operation. A conversion is considered valid
    when it produces a non-missing value.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    # Dunder methods
    # ==================================================================================
    def __str__(self) -> str:
        """Return the validator name formatted as a zero-argument call."""
        return super().__str__() + "()"

    # ==================================================================================
    def __call__(self, data: Series) -> tuple[Series, Series]:
        """Convert ``data`` and return its validity mask and converted values."""
        if data.empty:
            return self._return_empty(self._output_dtype, data.index)

        # ------------------------------------------------------------------------------
        ignored = self._ignored_mask(data)
        converted = self._convert(data)
        converted = converted.where(~ignored, self._fill_value)
        valid_mask = ignored | converted.notna()

        return self._return_validation_result(
            valid_mask, converted, error_message=f"Invalid data found in {self}:"
        )

    # ==================================================================================
    # Abstract properties and methods
    # ==================================================================================
    @property
    @abstractmethod
    def _output_dtype(self) -> type | str:
        """Dtype used when creating the converted Series for empty input."""
        pass

    # ==================================================================================
    @property
    @abstractmethod
    def _fill_value(self) -> Any:
        """Value assigned to entries excluded by the ignore rules."""
        pass

    # ==================================================================================
    @abstractmethod
    def _convert(self, data: Series) -> Series:
        """Convert ``data`` and represent unsuccessful conversions as missing."""
        pass
