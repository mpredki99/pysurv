# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from abc import ABC
from numbers import Number
from typing import Any, Callable, Iterable

import pandas as pd

from .constants import COMMENT
from .pysurv_validator import PySurvValidator
from .utils import PySurvValidatorMode


class ComparisonValidator(PySurvValidator, ABC):
    _operator: Callable[[Any, Any], Any]  # Operator defined in concrete child class

    def __init__(
        self,
        threshold: float,
        *,
        ignore: Iterable[Any] | Any = (pd.NA, COMMENT),
        mode: PySurvValidatorMode | str = PySurvValidatorMode.RAISE,
    ) -> None:
        super().__init__(ignore=ignore, mode=mode)
        self.threshold = threshold

    def __call__(self, data: pd.Series) -> tuple[pd.Series, pd.Series]:
        if data.empty:
            return self._return_empty(int, data.index)

        comparison = pd.Series(False, index=data.index)
        ignored = self._ignored_mask(data)

        try:
            comparison[~ignored] = self._operator(data[~ignored], self.threshold)
        except TypeError:
            return self._on_type_error(data)

        valid_mask = ignored | comparison

        return self._return_validation_result(valid_mask, data)

    def _on_type_error(self, data: pd.Series) -> tuple[pd.Series, pd.Series]:
        """Raise ValidationError or return result on TypeError during comparison."""
        mask = data.apply(lambda x: isinstance(x, Number))
        return self._return_validation_result(
            mask,
            data,
            error_message=f"Non-numeric data found in {self}:",
        )

    @property
    def threshold(self) -> float:
        return self._threshold

    @threshold.setter
    def threshold(self, value: float) -> None:
        self._threshold = float(value)
