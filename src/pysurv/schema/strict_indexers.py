# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from abc import ABC, abstractmethod
from typing import Any, Tuple

import pandas as pd


class StrictIndexer(ABC):

    _accessor_name: str  # Pandas accessor: 'loc', 'iloc', 'at', 'iat'

    def __init__(self, schema: "StrictSchema") -> None:
        self._schema = schema
        self._model = self._schema._model

    # ----------------------------------------------------------------------------------
    # Dunder methods
    # ----------------------------------------------------------------------------------
    def __getitem__(self, key: Any) -> Any:
        accessor = getattr(self._model, self._accessor_name)
        return accessor[key]

    def __setitem__(self, key: Any, value: Any) -> None:
        row_labels, col_labels = self._process_key(key)
        self._schema._assert_structure(set(row_labels), self._model.index, axis="row")
        self._schema._assert_structure(
            set(col_labels), self._model.columns, axis="column"
        )
        self._schema._assert_assignment(col_labels, value)
        self._set_value(key, value)

    def __delitem__(self, key: Any) -> None:
        raise KeyError(
            f"{self._schema.__class__.__name__} does not allow deleting item"
        )

    # ----------------------------------------------------------------------------------
    # Internal helpers
    # ----------------------------------------------------------------------------------
    def _get_accessor(self) -> Any:
        """Return the configured pandas indexer."""
        accessor = getattr(self._model, self._accessor_name, None)
        if accessor is None:
            raise AttributeError(f"Invalid accessor: {self._accessor_name}")
        return accessor

    def _set_value(self, key: Any, value: Any) -> None:
        """Set value on the model DataFrame using the proper accessor."""
        accessor = self._get_accessor()
        accessor[key] = value

    @abstractmethod
    def _process_key(self, key: Any) -> Tuple[list, list]:
        """Retrieve model row lables and column labels from key."""
        pass


# --------------------------------------------------------------------------------------
# Strict Block Indexer
# --------------------------------------------------------------------------------------
class StrictBlockIndexer(StrictIndexer):
    """Base class for strict loc/iloc indexers"""

    def _labels_from_result(self, tmp: pd.Series | pd.DataFrame) -> Tuple[set, set]:
        """Derive row and column label sets from an indexer result."""
        if isinstance(tmp, pd.Series):
            row_labels = (
                [tmp.name] if tmp.name in self._model.index else list(tmp.index)
            )
            col_labels = (
                [tmp.name] if tmp.name in self._model.columns else list(tmp.index)
            )
        else:
            row_labels = list(tmp.index)
            col_labels = list(tmp.columns)

        return row_labels, col_labels


# --------------------------------------------------------------------------------------
# Strict Loc Indexer
# --------------------------------------------------------------------------------------
class StrictLocIndexer(StrictBlockIndexer):
    _accessor_name = "loc"

    def _process_key(self, key: Any) -> Tuple[list, list]:
        tmp = self._get_accessor()[key]

        if isinstance(tmp, (pd.Series, pd.DataFrame)):
            return self._labels_from_result(tmp)

        row_labels = [key[0]]
        col_labels = [key[1]]
        return row_labels, col_labels


# --------------------------------------------------------------------------------------
# Strict Iloc Indexer
# --------------------------------------------------------------------------------------
class StrictIlocIndexer(StrictBlockIndexer):
    _accessor_name = "iloc"

    def _process_key(self, key: Any) -> Tuple[list, list]:
        tmp = self._get_accessor()[key]
        if isinstance(tmp, (pd.Series, pd.DataFrame)):
            return self._labels_from_result(tmp)

        row_labels = [self._model.index[key[0]]]
        col_labels = [self._model.columns[key[1]]]
        return row_labels, col_labels


# --------------------------------------------------------------------------------------
# Strict At Indexer
# --------------------------------------------------------------------------------------
class StrictAtIndexer(StrictIndexer):
    _accessor_name = "at"

    def _process_key(self, key: Tuple[Any, Any]) -> Tuple[list, list]:
        row_labels = [key[0]]
        col_labels = [key[1]]

        return row_labels, col_labels


# --------------------------------------------------------------------------------------
# Strict Iat Indexer
# --------------------------------------------------------------------------------------
class StrictIatIndexer(StrictIndexer):
    _accessor_name = "iat"

    def _process_key(self, key: Tuple[int, int]) -> Tuple[list, list]:
        row_labels = [self._model.index[key[0]]]
        col_labels = [self._model.columns[key[1]]]
        return row_labels, col_labels
