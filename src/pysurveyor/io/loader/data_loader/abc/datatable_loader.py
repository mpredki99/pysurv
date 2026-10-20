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
from pysurveyor.io.source.file_sources import CSVSource
from pysurveyor.typing.vars import DataTableType, SchemaType


# ======================================================================================
# Abstract base class for data-table loaders
# ======================================================================================
class DataTableLoader(TableLoader[DataTableType], Generic[DataTableType, SchemaType]):
    """Base loader for PySurveyor data tables.

    ``DataTableType`` identifies the concrete table returned by the loader.
    ``SchemaType`` identifies the schema accepted by the loader and forwarded
    to :meth:`from_csv`.

    The source dispatch implemented here supports :class:`CSVSource`.
    Concrete subclasses implement the actual CSV reading and construction of
    the data-table object.
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
        """Load a data table with a schema.

        This override only specializes the public signature for data-table
        loaders. Dispatch, source resolution, and forwarding remain owned by
        :meth:`TableLoader.__call__`.

        Parameters
        ----------
        source : Any
            Data source to load, resolved automatically when necessary.
        schema : SchemaType, optional
            Schema used to interpret or validate the data.
        config : IOConfig, optional
            Configuration used when resolving a path-like source. Concrete
            CSV loaders expect a ``CSVConfig`` instance.
        **kwargs
            Additional arguments forwarded to the concrete loader.

        Returns
        -------
        DataTableType
            Data-table object created by the concrete loader.
        """
        return super().__call__(source, schema=schema, config=config, **kwargs)

    # ==================================================================================
    def from_source(
        self,
        source: Source,
        schema: SchemaType | None = None,
        **kwargs: Any,
    ) -> DataTableType:
        """Load a data table from a supported resolved source.

        Parameters
        ----------
        source : Source
            Resolved source describing the input data.
        schema : SchemaType, optional
            Schema used to interpret or validate the data.
        **kwargs
            Additional arguments forwarded to :meth:`from_csv`.

        Returns
        -------
        DataTableType
            Data-table object created by the concrete loader.

        Raises
        ------
        TypeError
            If ``source`` is not a supported source type.
        """
        if isinstance(source, CSVSource):
            return self.from_csv(
                source.path,
                schema=schema,
                config=source.config,
                **kwargs,
            )

        raise TypeError(f"Unsupported source type: {type(source)}")

    # ==================================================================================
    @abstractmethod
    def from_csv(
        self,
        path: PathLike[str] | str,
        schema: SchemaType | None = None,
        config: CSVConfig | None = None,
        **kwargs: Any,
    ) -> DataTableType:
        """Load a data table from a CSV file."""
        pass
