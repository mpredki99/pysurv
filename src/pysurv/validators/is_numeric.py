# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import numpy as np
import pandas as pd

from .pysurv_validator import PySurvValidator


class IsNumeric(PySurvValidator):
    def __call__(self, data: pd.Series) -> tuple[pd.Series, pd.Series]:
        """Validate numeric values in the Series. Normalize ignored values."""
        if data.empty:
            return self.return_empty(float, data.index)

        ignored = self.ignored_mask(data)

        converted = pd.to_numeric(data, errors="coerce")
        converted = converted.where(~ignored, np.nan)

        valid_mask = ignored | converted.notna()

        return self.return_validation_result(
            valid_mask, converted, error_message=f"Invalid data found in {self}:"
        )
