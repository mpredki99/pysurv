# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from typing import Any

from pysurveyor.backend.backend import DataFrame, Index
from pysurveyor.io.abc.pysurveyor_io import PySurveyorIO
from pysurveyor.typing.vars import LoaderType


# ======================================================================================
# Abstract base class for DataFrame-like PySurveyor objects
# ======================================================================================
class Table(PySurveyorIO[LoaderType]):
    """Base class for PySurveyor objects backed by tabular data.

    ``Table`` wraps a pandas :class:`DataFrame` and exposes a minimal, explicit
    interface for accessing and modifying its underlying data.

    The public interface includes item access and assignment, standard string
    representations, index and column labels, and indexing accessors.  Additional
    DataFrame functionality may be exposed explicitly as the public API evolves.

    Parameters
    ----------
    data : DataFrame
        Tabular data represented by the object.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    # Dunder methods
    # ==================================================================================
    def __init__(self, data: DataFrame) -> None:
        """Initialize a table with its underlying data."""
        super().__init__()
        self._data = data

    # ==================================================================================
    def __getitem__(self, key: Any) -> Any:
        """Return an item from the underlying DataFrame."""
        return self._data[key]

    # ==================================================================================
    def __setitem__(self, key: Any, value: Any) -> None:
        """Set an item in the underlying DataFrame."""
        self._data[key] = value

    # ==================================================================================
    def __repr__(self) -> str:
        """Return the representation of the underlying data."""
        return repr(self._data)

    # ==================================================================================
    def __str__(self) -> str:
        """Return the string representation of the underlying data."""
        return str(self._data)

    # ==================================================================================
    # Properties
    # ==================================================================================
    @property
    def index(self) -> Index:
        """Return the index of the underlying data."""
        return self._data.index

    # ==================================================================================
    @property
    def columns(self) -> Index:
        """Return the column names of the table."""
        return self._data.columns

    # ==================================================================================
    # Indexers
    # ==================================================================================
    @property
    def loc(self):
        """Provide label-based indexing for rows and columns."""
        return self._data.loc

    # ==================================================================================
    @property
    def iloc(self):
        """Provide integer-based indexing for rows and columns."""
        return self._data.iloc

    # ==================================================================================
    @property
    def at(self):
        """Provide label-based scalar access for rows and columns."""
        return self._data.at

    # ==================================================================================
    @property
    def iat(self):
        """Provide integer-based scalar access for rows and columns."""
        return self._data.iat
