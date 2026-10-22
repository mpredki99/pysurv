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
        self._inclusive = self._validate_inclusive(inclusive)

        # ------------------------------------------------------------------------------
        if self._minimum > self._maximum:
            raise ValueError("minimum cannot be greater than maximum")

        # ------------------------------------------------------------------------------
        self._validator = self._build_validator()

    # ==================================================================================
    def __str__(self) -> str:
        """Return the validator name, limits, and inclusivity mode."""
        return (
            super().__str__()
            + f"({self.minimum}, {self.maximum}, inclusive={self.inclusive})"
        )

    # ==================================================================================
    def __call__(self, data: Series) -> tuple[Series, Series]:
        """Validate ``data`` using the configured numeric range."""
        return self._validator(data)

    # ==================================================================================
    # Properties
    # ==================================================================================
    @Validator.mode.setter
    def mode(self, value: ValidatorMode | str) -> None:
        """Set this validator's mode and propagate it to its range checks."""
        setter = Validator.mode.fset

        # ------------------------------------------------------------------------------
        if setter is None:
            raise AttributeError("Validator.mode does not define a setter")
        setter(self, value)

        # ------------------------------------------------------------------------------
        if hasattr(self, "_validator"):
            self._validator.mode = self._mode

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
            self._refresh_validator()

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
            self._refresh_validator()

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
            self._refresh_validator()

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
    def _build_validator(self) -> Validator:
        """Build the range validator from the current configuration."""
        left_cls, right_cls = _RANGE_OPERATORS[self._inclusive]

        # ------------------------------------------------------------------------------
        validator = left_cls(self._minimum, ignore=self.ignore) & right_cls(
            self._maximum, ignore=self.ignore
        )
        validator.mode = self._mode

        # ------------------------------------------------------------------------------
        return validator

    # ==================================================================================
    def _refresh_validator(self) -> None:
        """Rebuild the range validator after a configuration change."""
        if hasattr(self, "_validator"):
            self._validator = self._build_validator()
