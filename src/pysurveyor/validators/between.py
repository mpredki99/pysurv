# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from collections.abc import Iterable
from typing import Any, Literal, cast

from pysurveyor.backend.backend import NA, Series
from pysurveyor.validators.abc.validator import Validator
from pysurveyor.validators.comparison_validators import (
    Greater,
    GreaterOrEqual,
    Less,
    LessOrEqual,
)
from pysurveyor.validators.conversion_validators import IsNumeric
from pysurveyor.validators.utils.constants import COMMENT
from pysurveyor.validators.utils.validator_mode import ValidatorMode

# ======================================================================================
Inclusivity = Literal["both", "left", "right", "neither"]


# ======================================================================================
_RANGE_OPERATORS = {
    "both": (GreaterOrEqual, LessOrEqual),
    "left": (GreaterOrEqual, Less),
    "right": (Greater, LessOrEqual),
    "neither": (Greater, Less),
}


# ======================================================================================
# Between validator
# ======================================================================================
class Between(Validator):
    """Validate values against an inclusive or exclusive numeric range."""

    # ==================================================================================
    # Dunder methods
    # ==================================================================================
    def __init__(
        self,
        minimum: float,
        maximum: float,
        *,
        inclusive: Inclusivity = "left",
        ignore: Iterable[Any] | Any = (NA, COMMENT),
        mode: ValidatorMode | str = ValidatorMode.RAISE,
    ) -> None:
        """Initialize the range limits and validation configuration."""
        super().__init__(ignore=ignore, mode=mode)

        # ------------------------------------------------------------------------------
        self._minimum = float(minimum)
        self._maximum = float(maximum)

        # ------------------------------------------------------------------------------
        if self._minimum > self._maximum:
            raise ValueError("minimum cannot be greater than maximum")

        # ------------------------------------------------------------------------------
        self._inclusive = self._validate_inclusive(inclusive)
        self._is_numeric = IsNumeric(mode="return")

    def __str__(self) -> str:
        """Return the validator name, limits, and inclusivity mode."""
        return (
            super().__str__()
            + f"({self.minimum}, {self.maximum}, inclusive={self.inclusive})"
        )

    # ==================================================================================
    def __call__(self, data: Series) -> tuple[Series, Series]:
        """Validate ``data`` using the configured numeric range."""
        if data.empty:
            return self._return_empty(int, data.index)

        # ------------------------------------------------------------------------------
        ignored = self._ignored_mask(data)
        numeric = self._is_numeric(data)[0]
        comparable = ~ignored & numeric

        # ------------------------------------------------------------------------------
        valid_mask = ignored.copy()

        if comparable.any():
            valid_mask.loc[comparable] = self._validator(data.loc[comparable])[0]

        # ------------------------------------------------------------------------------
        error_message = (
            f"Non-numeric data found:"
            if (~ignored & ~numeric).any()
            else f"Data out of range in:"
        )

        # ------------------------------------------------------------------------------
        return self._return_validation_result(
            valid_mask,
            data,
            error_message=error_message,
        )

    # ==================================================================================
    # Properties
    # ==================================================================================
    @property
    def inclusive(self) -> Inclusivity:
        """Return the range boundary inclusivity mode."""
        return cast(Inclusivity, self._inclusive)

    @inclusive.setter
    def inclusive(self, value: Inclusivity) -> None:
        """Set the range boundary inclusivity mode."""
        value = self._validate_inclusive(value)

        # ------------------------------------------------------------------------------
        if value != self._inclusive:
            self._inclusive = value

    # ==================================================================================
    @property
    def minimum(self) -> float:
        """Return the lower range limit."""
        return self._minimum

    @minimum.setter
    def minimum(self, value: float) -> None:
        """Set the lower range limit as a float."""
        value = float(value)

        # ------------------------------------------------------------------------------
        if value > self._maximum:
            raise ValueError("minimum cannot be greater than maximum")

        # ------------------------------------------------------------------------------
        if value != self._minimum:
            self._minimum = value

    # ==================================================================================
    @property
    def maximum(self) -> float:
        """Return the upper range limit."""
        return self._maximum

    @maximum.setter
    def maximum(self, value: float) -> None:
        """Set the upper range limit as a float."""
        value = float(value)

        # ------------------------------------------------------------------------------
        if value < self._minimum:
            raise ValueError("minimum cannot be greater than maximum")

        # ------------------------------------------------------------------------------
        if value != self._maximum:
            self._maximum = value

    # ==================================================================================
    # Private methods
    # ==================================================================================
    @staticmethod
    def _validate_inclusive(value: str) -> Inclusivity:
        """Validate the range boundary inclusivity mode."""
        if value not in _RANGE_OPERATORS:
            raise ValueError(
                "inclusive must be one of: 'both', 'left', 'right', 'neither'"
            )
        return cast(Inclusivity, value)

    # ==================================================================================
    @property
    def _validator(self) -> Validator:
        """Return a composite validator based on the configured range and inclusivity."""
        left_op, right_op = _RANGE_OPERATORS[self._inclusive]

        # ------------------------------------------------------------------------------
        left = left_op(self._minimum, ignore=self._ignore)
        right = right_op(self._maximum, ignore=self._ignore)

        # ------------------------------------------------------------------------------
        validator = left & right
        validator.mode = ValidatorMode.RETURN

        return validator
