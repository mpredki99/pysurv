# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from abc import abstractmethod
from os import PathLike
from typing import Any, Generic

from pysurveyor.io.config.abc.config import IOConfig
from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.loader.abc.table_loader import TableLoader
from pysurveyor.io.source.abc.source import Source
from pysurveyor.typing.vars import DataTableType, SchemaType


# ======================================================================================
# Abstract base class for data-table loaders
# ======================================================================================
class DataTableLoader(TableLoader[DataTableType], Generic[DataTableType, SchemaType]):
    """Base loader for PySurveyor data-table objects.

    This abstract base class defines the interface for loaders that construct
    data tables associated with schemas. It resolves supported source types
    and delegates format-specific loading to concrete implementations.

    The generic parameters specify the type of table returned by the loader
    and the schema type accepted by its loading methods.

    Parameters
    ----------
    DataTableType
        Concrete data-table type returned by the loader.
    SchemaType
        Schema type accepted by the loader and passed to the format-specific
        loading implementation.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    def __call__(
        self,
        source: Any,
        *,
        schema: SchemaType | None = None,
        config: IOConfig | None = None,
        **kwargs: Any,
    ) -> DataTableType:
        """Load a data table from a source.

        Delegates source resolution and loading orchestration to
        :meth:`TableLoader.__call__`, forwarding the schema, configuration,
        and additional keyword arguments.

        Parameters
        ----------
        source : Any
            Data source to load from, resolved automatically.
        schema : SchemaType or None, default=None
            Schema to associate with or use when interpreting the loaded data.
            If ``None``, schema selection is delegated to the concrete loader
            or data-table implementation.
        config : IOConfig, optional
            Configuration of the reading the source. If ``source`` is already
            a :class:`Source`, its config object is used.
        **kwargs : Any
            Additional keyword arguments forwarded to the inherited loading
            implementation.

        Returns
        -------
        DataTableType
            Data-table instance created by the concrete loader.
        """
        return super().__call__(source, schema=schema, config=config, **kwargs)

    # ==================================================================================
    def from_source(
        self,
        source: Source,
        *,
        schema: SchemaType | None = None,
        **kwargs: Any,
    ) -> DataTableType:
        """Load a data table from a supported source.

        Parameters
        ----------
        source : Source
            Source describing the input data and its configuration.
        schema : SchemaType, optional
            Schema used to interpret or validate the data.
        **kwargs
            Additional arguments forwarded to source specific loading methods.

        Returns
        -------
        DataTableType
            Data-table object created by the concrete loader.

        Raises
        ------
        TypeError
            If ``source`` is not a supported source type.
        """
        return super().from_source(source, schema=schema, **kwargs)
