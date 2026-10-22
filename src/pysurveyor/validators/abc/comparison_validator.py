# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from collections.abc import Iterable
from numbers import Number
from typing import Any, Callable, cast

from pysurveyor.abc.pysurveyor_abc import ABSTRACT_CLASS_ATTRIBUTE
from pysurveyor.backend.backend import NA, Series
from pysurveyor.validators.abc.validator import Validator
from pysurveyor.validators.utils.constants import COMMENT
from pysurveyor.validators.utils.validator_mode import ValidatorMode


# ======================================================================================
# Comparison validator base class
# ======================================================================================
class ComparisonValidator(Validator):
    """Base class for validators that compare numeric values with a threshold.

    Concrete validators provide the comparison operator. Values matching the
    configured ignore rules are considered valid without being compared.
    """

    __pysurveyor_abstract__ = True

    _operator = cast(Callable, ABSTRACT_CLASS_ATTRIBUTE)
    _operator_error_str = cast(str, ABSTRACT_CLASS_ATTRIBUTE)

    # ==================================================================================
    # Dunder methods
    # ==================================================================================
    def __init__(
        self,
        threshold: float,
        *,
        ignore: Iterable[Any] | Any = (NA, COMMENT),
        mode: ValidatorMode | str = ValidatorMode.RAISE,
    ) -> None:
        """Initialize the comparison threshold and validation configuration."""
        super().__init__(ignore=ignore, mode=mode)
        self.threshold = threshold

    # ==================================================================================
    def __str__(self) -> str:
        """Return the validator name and configured threshold."""
        return super().__str__() + f"({self.threshold})"

    # ==================================================================================
    def __call__(self, data: Series) -> tuple[Series, Series]:
        """Compare ``data`` with the threshold and return validity results.

        Ignored values are marked valid. If comparison raises ``TypeError``,
        non-numeric values are reported as invalid through ``_on_type_error``.
        """
        if data.empty:
            return self._return_empty(int, data.index)

        # ------------------------------------------------------------------------------
        comparison = Series(False, index=data.index)
        ignored = self._ignored_mask(data)

        try:
            comparison[~ignored] = self._operator(data[~ignored], self.threshold)
        except TypeError:
            return self._on_type_error(data)

        valid_mask = ignored | comparison

        return self._return_validation_result(
            valid_mask,
            data,
            error_message=f"Found values {self._operator_error_str} {self.threshold}:",
        )

    # ==================================================================================
    # Properties
    # ==================================================================================
    @property
    def threshold(self) -> float:
        """Return the numeric comparison threshold."""
        return self._threshold

    @threshold.setter
    def threshold(self, value: float) -> None:
        """Set the comparison threshold as a float."""
        self._threshold = float(value)

    # ==================================================================================
    # Private methods
    # ==================================================================================
    def _on_type_error(self, data: Series) -> tuple[Series, Series]:
        """Handle non-numeric values after a failed comparison."""
        mask = self._ignored_mask(data) | data.apply(lambda x: isinstance(x, Number))
        return self._return_validation_result(
            mask,
            data,
            error_message=f"Non-numeric data found in {self}:",
        )
