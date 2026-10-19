# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from pysurveyor.backend.backend import DataFrame
from pysurveyor.io.abc.pysurveyor_io import PySurveyorIO
from pysurveyor.typing.vars import LoaderType


# ======================================================================================
# Abstract base class for DataFrame-like PySurveyor objects
# ======================================================================================
class Table(PySurveyorIO[LoaderType]):
    """Base class for PySurveyor objects backed by tabular data.

    ``Table`` stores the underlying :class:`DataFrame` and provides its
    standard string representations.

    Parameters
    ----------
    data : DataFrame
        Tabular data represented by the object.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    def __init__(self, data: DataFrame) -> None:
        """Initialize a table with its underlying data."""
        super().__init__()
        self._data = data

    # ==================================================================================
    def __repr__(self) -> str:
        """Return the representation of the underlying data."""
        return repr(self._data)

    # ==================================================================================
    def __str__(self) -> str:
        """Return the string representation of the underlying data."""
        return str(self._data)
