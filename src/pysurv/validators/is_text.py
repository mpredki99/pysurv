# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pandas as pd

from .conversion_validator import ConversionValidator


class IsText(ConversionValidator):
    @property
    def _output_dtype(self) -> str:
        return "str"

    @property
    def _fill_value(self) -> pd.NA:
        return pd.NA

    def _convert(self, data: pd.Series) -> pd.Series:
        return data.astype("str")
