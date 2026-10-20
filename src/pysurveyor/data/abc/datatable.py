# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from typing import Any, Generic, Self, cast

from pysurveyor.abc.pysurveyor_abc import ABSTRACT_CLASS_ATTRIBUTE
from pysurveyor.backend.backend import DataFrame
from pysurveyor.io.abc.table import Table
from pysurveyor.typing.vars import LoaderType, SchemaType


# ======================================================================================
# Abstract base class data tables with associated schemas
# ======================================================================================
class DataTable(Table[LoaderType], Generic[SchemaType, LoaderType]):
    """Base class for data tables associated with a typed schema.

    The schema describes the structure and metadata of the table's data.
    Concrete subclasses specify the schema type through ``SchemaType`` and
    the loader type exposed by the inherited ``load`` accessor through
    ``LoaderType``.

    Parameters
    ----------
    data : DataFrame
        Tabular data represented by the object.
    schema : SchemaType, optional
        Schema associated with the table. If ``None``, a deep copy of the
        concrete class's ``default_schema`` is used.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    default_schema = cast(SchemaType, ABSTRACT_CLASS_ATTRIBUTE)

    # ==================================================================================
    def __init__(self, data: DataFrame, schema: SchemaType | None = None) -> None:
        """Initialize a data table and associate its schema."""
        super().__init__(data)
        self._schema = schema or self.__class__.default_schema.copy(deep=True)

    # ==================================================================================
    @property
    def schema(self) -> SchemaType:
        """Return the schema associated with this data table."""
        return self._schema

    # ==================================================================================
    def copy(self, deep: bool = True) -> Self:
        """Return a copy of the table, its schema, and copy options.

        Parameters
        ----------
        deep : bool, default=True
            Whether to copy the table data and schema deeply.

        Returns
        -------
        Self
            A copy of the same concrete data-table type.
        """
        return self.__class__(
            self._data.copy(deep=deep),
            self._schema.copy(deep=deep),
            **self._copy_kwargs(),
        )

    # ==================================================================================
    def _copy_kwargs(self) -> dict[str, Any]:
        """Return additional keyword arguments for constructing a table copy.

        Concrete subclasses can override this method to provide constructor
        arguments that must be preserved when copying an instance.

        Returns
        -------
        dict[str, Any]
            Additional keyword arguments passed to the concrete class's constructor.
        """
        return {}
