# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from typing import Any, Iterable

import pandas as pd

from .pysurv_validator import PySurvValidator, PySurvValidatorMode


class GreaterOrEqual(PySurvValidator):
    def __init__(
        self,
        threshold: float = 0.0,
        *,
        ignore: Iterable[Any] | Any = pd.NA,
        mode: PySurvValidatorMode | str = PySurvValidatorMode.RAISE,
    ) -> None:
        """Configure ignore values set and threshold value."""
        super().__init__(ignore=ignore, mode=mode)
        self.threshold = threshold

    def __str__(self) -> str:
        return f"GreaterOrEqual {self.threshold} (ignore={self.ignore!r})"

    def __call__(self, data: pd.Series) -> tuple[pd.Series, pd.Series]:
        if data.empty:
            return self.return_empty(float, data.index)

        ignored = self.ignored_mask(data)

        try:
            comparison = data >= self.threshold
        except TypeError:
            return self.return_validation_result(
                data.apply(lambda x: isinstance(x, (float, int))),
                data,
                error_message=f"Non-numeric data found in {self}:",
            )

        valid_mask = ignored | comparison

        return self.return_validation_result(
            valid_mask, data, error_message=f"Invalid data found in {self}:"
        )

    @property
    def threshold(self) -> float:
        """Return threshold value."""
        return self._threshold

    @threshold.setter
    def threshold(self, value: float) -> None:
        """Set threshold value and ensure numeric type."""
        try:
            self._threshold = float(value)
        except ValueError:
            raise ValueError(
                f"{self.__class__.__name__}.threshold must be numeric value."
            )
