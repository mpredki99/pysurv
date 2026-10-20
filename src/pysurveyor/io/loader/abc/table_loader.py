# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from typing import Any

from pysurveyor.io.config.abc.config import IOConfig
from pysurveyor.io.loader.abc.constructor_bound_loader import ConstructorBoundLoader
from pysurveyor.io.source.abc.source import Source
from pysurveyor.io.source.file_sources import CSVSource
from pysurveyor.io.source.resolver import resolve_source
from pysurveyor.typing.vars import DomainClassType


# ======================================================================================
# All Table objects loader
# ======================================================================================
class TableLoader(ConstructorBoundLoader[DomainClassType]):
    """Base loader for PySurveyor ``Table`` objects.

    ``TableLoader`` extends :class:`ConstructorBoundLoader`, so each concrete
    loader is bound to the ``Table`` class constructor.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    def __call__(
        self,
        source: Any,
        *,
        config: IOConfig | None = None,
        **kwargs: Any,
    ) -> DomainClassType:
        """Load a domain object from a source or path-like input.

        Parameters
        ----------
        source : Any
            Data source to load from, resolved automatically.
        config : IOConfig, optional
            Configuration of the reading the source. If ``source`` is already
            a :class:`Source`, its config object is used.
        **kwargs
            Arguments forwarded to subclass implementations.

        Returns
        -------
        DomainClassType
            Domain object created by the concrete loader.
        """
        if isinstance(source, Source):
            resolved_source = source
        else:
            resolved_source = resolve_source(source, config=config)

        return self.from_source(resolved_source, **kwargs)

    # ==================================================================================
    def from_source(
        self,
        source: Source,
        **kwargs: Any,
    ) -> DomainClassType:
        """Load a domain object from a supported source.

        Parameters
        ----------
        source : Source
            Source describing the domain object input and its configuration.
        **kwargs
            Additional arguments forwarded to source specific loading methods.

        Returns
        -------
        DomainClassType
            Domain object created by the concrete loader.

        Raises
        ------
        TypeError
            If ``source`` is not a supported source type.
        """
        if isinstance(source, CSVSource):
            return self.from_csv(source.path, config=source.config, **kwargs)

        raise TypeError(f"Unsupported source type: {type(source)}")
