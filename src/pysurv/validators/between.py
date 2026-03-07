# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from typing import Any, Iterable

import pandas as pd

from .constants import COMMENT
from .greater import Greater
from .greater_equal import GreaterOrEqual
from .less import Less
from .less_equal import LessOrEqual
from .pysurv_validator import PySurvValidator
from .utils import PySurvValidatorMode


class Between(PySurvValidator):
    def __init__(
        self,
        minimum: float,
        maximum: float,
        *,
        inclusive: str = "left",
        ignore: Iterable[Any] | Any = (pd.NA, COMMENT),
        mode: PySurvValidatorMode | str = PySurvValidatorMode.RAISE,
    ) -> None:
        super().__init__(ignore=ignore, mode=mode)

        self.minimum = minimum
        self.maximum = maximum
        self.inclusive = inclusive

        if inclusive == "both":
            left = GreaterOrEqual(self.minimum, ignore=ignore)
            right = LessOrEqual(self.maximum, ignore=ignore)

        elif inclusive == "left":
            left = GreaterOrEqual(self.minimum, ignore=ignore)
            right = Less(self.maximum, ignore=ignore)

        elif inclusive == "right":
            left = Greater(self.minimum, ignore=ignore)
            right = LessOrEqual(self.maximum, ignore=ignore)

        else:
            left = Greater(self.minimum, ignore=ignore)
            right = Less(self.maximum, ignore=ignore)

        self._validator = left & right
        self._validator.mode = self.mode

    def __call__(self, data: pd.Series) -> tuple[pd.Series, pd.Series]:
        return self._validator(data)

    @property
    def inclusive(self) -> str:
        return self._inclusive

    @inclusive.setter
    def inclusive(self, value: str) -> None:
        if value not in {"both", "left", "right", "neither"}:
            raise ValueError(
                "inclusive must be one of: 'both', 'left', 'right', 'neither'"
            )
        self._inclusive = value

    @property
    def minimum(self) -> float:
        return self._minimum

    @minimum.setter
    def minimum(self, value: float) -> None:
        self._minimum = float(value)

    @property
    def maximum(self) -> float:
        return self._maximum

    @maximum.setter
    def maximum(self, value: float) -> None:
        value = float(value)
        if value < self._minimum:
            raise ValueError("minimum cannot be greater than maximum")
        self._maximum = value
