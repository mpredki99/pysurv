# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from collections.abc import Iterable
from itertools import cycle
from typing import Any

import pandas as pd

from pysurv.typing.typing import PySurvUnit
from pysurv.validators.pysurv_validator import PySurvValidator

from .field_compilers import compile_unit, compile_validator
from .pysurv_schema import PySurvSchema
from .strict_indexers import (
    StrictAtIndexer,
    StrictIatIndexer,
    StrictIlocIndexer,
    StrictLocIndexer,
)
from .utils import SchemaValidationMode


class StrictSchema(PySurvSchema):
    """
    PySurvSchema subclass that enforce its structure, so that no schema model fields
    can be added, removed or renamed. `Required`, `Unit`, `Validator` and `Description`
    values still can be changed with strict type checking.
    """

    _column_types = {
        "required": bool,
        "unit": PySurvUnit,
        "validator": PySurvValidator,
        "description": str,
    }

    _indexers = {
        "loc": StrictLocIndexer,
        "iloc": StrictIlocIndexer,
        "at": StrictAtIndexer,
        "iat": StrictIatIndexer,
    }

    _forbidden_methods = {
        "where",
        "mask",
        "assign",
        "update",
        "insert",
        "pop",
        "drop",
        "rename",
        "set_index",
        "reset_index",
    }

    def __init__(
        self,
        model: pd.DataFrame,
        validation_mode: SchemaValidationMode | str = SchemaValidationMode.LAZY,
    ) -> None:
        object.__setattr__(self, "_initialized", False)
        super().__init__(model, validation_mode)
        self._row_set = set(self._model.index)
        self._col_set = set(self._model.columns)
        object.__setattr__(self, "_initialized", True)
        self._assert_init_type()

    # ----------------------------------------------------------------------------------
    # Dunder methods
    # ----------------------------------------------------------------------------------
    def __getattr__(self, name: str) -> Any:
        if name in self._indexers:
            return self._indexers[name](self)

        attr = getattr(self._model, name, None)
        if name not in self._forbidden_methods and attr is not None:
            return attr

        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    def __setattr__(self, name: str, value: Any) -> None:
        if object.__getattribute__(self, "_initialized") and not hasattr(self, name):
            raise AttributeError(f"'{self.__class__.__name__}' cannot be modified")
        super().__setattr__(name, value)

    def __delattr__(self, name: str) -> None:
        raise AttributeError(
            f"'{type(self).__name__}' does not allow deleting attribute '{name}'"
        )

    def __getitem__(self, key: Any) -> Any:
        return self._model[key]

    def __setitem__(self, key: Any, value: Any) -> None:
        tmp = self._model[key]
        row_labels = list(tmp.index)
        col_labels = [tmp.name] if isinstance(tmp, pd.Series) else list(tmp.columns)

        self._assert_structure(set(row_labels), self._row_set, "row")
        self._assert_structure(set(col_labels), self._col_set, "column")
        self._assert_assignment(col_labels, value)

        self._model[key] = value

    # ----------------------------------------------------------------------------------
    # Ensure strict model's structure
    # ----------------------------------------------------------------------------------
    def _assert_structure(
        self, labels: set[str], existing: set[str], axis: str
    ) -> None:
        """Prevent from creating new rows or columns."""
        if not labels.issubset(existing):
            new = [label for label in labels if label not in existing]
            raise KeyError(f"Cannot create new {axis}: {", ".join(new)}")

    def _assert_assignment(self, col_labels: list[str], value: Any) -> None:
        """Validate assigned values are correct type and prevent col number mismatch."""
        is_scalar = not isinstance(value, Iterable) or isinstance(value, (str, bytes))
        values = [value] if is_scalar else list(value)

        self._assert_column_mismatch(col_labels, values)
        self._assert_types(col_labels, values)

    def _assert_column_mismatch(self, col_labels: list[str], values: list[Any]) -> None:
        """Ensure number of columns match number of values."""
        n_values = len(values)
        n_cols = len(col_labels)

        if n_values > 1 and n_values != n_cols:
            raise ValueError(
                f"Column number mismatch: {n_values} values for {n_cols} columns"
            )

    def _assert_types(self, col_labels: list[str], values: list[Any]) -> None:
        """Ensure number assigned values are correct type."""
        for col, val in zip(col_labels, cycle(values)):
            expected_type = self._column_types.get(col)

            if expected_type is None:
                raise ValueError(f"Unknown column '{col}' in schema")

            if not isinstance(val, expected_type):
                raise TypeError(
                    f"Column '{col}' expects '{expected_type.__name__}', got '{type(val).__name__}'"
                )

    def _assert_init_type(self):
        """Ensure the schema was created with proper initial column types."""
        for col in self._col_set:
            expected_type = self._column_types.get(col)

            invalid_values = ~self._model[col].apply(
                lambda x: True if pd.isna(x) else isinstance(x, expected_type)
            )
            if invalid_values.any():
                raise ValueError(
                    f"Wrong dtype in {col} detected: "
                    f"{list(self._model[col][invalid_values])}"
                )
