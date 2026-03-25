# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from collections import Counter
from pathlib import Path

from pysurv.schema.utils import get_models_dir

from .strict_schema import StrictSchema


class PySurvTableSchema(StrictSchema):
    # ----------------------------------------------------------------------------------
    # Validate given column subtypes attributes
    # ----------------------------------------------------------------------------------
    def _assert_column_subsets(self, *columns) -> None:
        """Validate provided column subtypes."""
        self._assert_column_subsets_in_schema(columns)
        self._assert_column_subsets_unique(columns)

    def _assert_column_subsets_in_schema(self, columns) -> None:
        """Ensure all provided column labels are in schema model."""
        schema_index = self._model.index
        for cols in columns:
            cols_set = set(cols)
            if not cols_set.issubset(schema_index):
                raise ValueError(
                    "Given columns are not in schema: "
                    f"{[col for col in cols if col not in schema_index]}"
                )

    def _assert_column_subsets_unique(self, columns) -> None:
        """Ensure all provided column labels are unique."""
        final_list = [col for cols in columns for col in cols]
        final_set = set(final_list)

        if len(final_list) != len(final_set):
            counts = Counter(final_list)
            duplicates = [col for col, count in counts.items() if count > 1]
            raise ValueError(f"Given column subsets are not unique: {duplicates}")

    # ----------------------------------------------------------------------------------
    # Internal helpers
    # ----------------------------------------------------------------------------------
    def _get_model_file_path(self, file_name: str) -> Path:
        """Return path to given file in models directory."""
        models_dir = get_models_dir()
        model_path = models_dir / file_name

        if not model_path.is_file():
            raise FileNotFoundError(
                f"{self.__class__.__name__} model file not found in: {model_path}"
            )

        return model_path
