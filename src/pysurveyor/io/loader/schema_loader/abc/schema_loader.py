# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from abc import abstractmethod
from os import PathLike
from typing import Any

from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.loader.abc.table_loader import TableLoader
from pysurveyor.io.source.abc.source import Source
from pysurveyor.io.source.file_sources import CSVSource
from pysurveyor.typing.vars import SchemaType


# ======================================================================================
# Abstract base class for schema loaders
# ======================================================================================
class SchemaLoader(TableLoader[SchemaType]):
    """Base loader for schema objects.

    This abstract base class defines the interface for loading schema objects and
    provides dispatch for supported source types. Source resolution and the public
    loading interface are inherited from :class:`TableLoader`.

    Type Parameters
    ---------------
    SchemaType
        Concrete schema type returned by the loader.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    def from_source(
        self,
        source: Source,
        **kwargs: Any,
    ) -> SchemaType:
        """Load a schema from a supported source.

        Parameters
        ----------
        source : Source
            Source describing the schema input and its configuration.
        **kwargs
            Additional arguments forwarded to source specific loading methods.

        Returns
        -------
        SchemaType
            Schema object created by the concrete loader.

        Raises
        ------
        TypeError
            If ``source`` is not a supported source type.
        """
        if isinstance(source, CSVSource):
            return self.from_csv(source.path, config=source.config, **kwargs)

        raise TypeError(f"Unsupported source type: {type(source)}")

    # ==================================================================================
    @abstractmethod
    def from_csv(
        self,
        path: PathLike[str] | str,
        config: CSVConfig | None = None,
        **kwargs: Any,
    ) -> SchemaType:
        """Load a schema object from a CSV file."""
        pass
