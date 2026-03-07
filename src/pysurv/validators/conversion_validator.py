# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from abc import ABC, abstractmethod
from typing import Any

import pandas as pd

from .pysurv_validator import PySurvValidator


class ConversionValidator(PySurvValidator, ABC):
    """Base for validators that convert values."""

    @property
    @abstractmethod
    def _output_dtype(self) -> type | str: ...

    @property
    @abstractmethod
    def _fill_value(self) -> Any:
        """Value to use for ignored entries."""
        ...

    @abstractmethod
    def _convert(self, data: pd.Series) -> pd.Series:
        """Convert data dtype."""
        ...

    def __call__(self, data: pd.Series) -> tuple[pd.Series, pd.Series]:
        if data.empty:
            return self._return_empty(self._output_dtype, data.index)

        ignored = self._ignored_mask(data)
        converted = self._convert(data)
        converted = converted.where(~ignored, self._fill_value)
        valid_mask = ignored | converted.notna()

        return self._return_validation_result(
            valid_mask, converted, error_message=f"Invalid data found in {self}:"
        )
