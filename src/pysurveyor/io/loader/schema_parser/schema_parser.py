# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from collections.abc import Iterable, Iterator
from os import PathLike
from pathlib import Path
import re
from typing import Any

from pysurveyor.backend.backend import read_csv
from pysurveyor.io.config.abc.config import IOConfig
from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.loader.abc.object_bound_loader import ObjectBoundLoader
from pysurveyor.io.source.abc.source import Source
from pysurveyor.io.source.file_sources import CSVSource
from pysurveyor.io.source.resolver import resolve_source
from pysurveyor.typing.vars import SchemaType


# ======================================================================================
# Schema parser
# ======================================================================================
class SchemaParser(ObjectBoundLoader[SchemaType]):

    # ==================================================================================
    def __call__(self, source: Any, *, config: IOConfig | None = None) -> SchemaType:
        """Parse a schema from a supported source.

        Parameters
        ----------
        source : Any
            Data source to load from, resolved automatically.
        config : IOConfig, optional
            Configuration of the reading the source. If ``source`` is already
            a :class:`Source`, its config object is used.

        Returns
        -------
        SchemaType
            Schema object with fields configured according to the header.
        """
        if isinstance(source, Source):
            resolved_source = source
        else:
            resolved_source = resolve_source(source, config=config)

        return self.from_source(resolved_source)

    # ==================================================================================
    def from_source(self, source: Source) -> SchemaType:
        """Dispatch a source to the schema parser."""
        if isinstance(source, CSVSource):
            return self.from_csv(source.path, config=source.config)

        raise TypeError(f"Unsupported source type: {type(source)}")

    # ==================================================================================
    def from_csv(
        self, path: PathLike[str] | str, config: CSVConfig | None = None
    ) -> SchemaType:
        """Read a CSV header and infer a schema from its annotations."""
        path = Path(path).expanduser().resolve()

        if not path.is_file():
            raise FileNotFoundError(f"CSV file not found: {path}")

        config = config or CSVConfig()

        read_kwargs = config.kwargs
        read_kwargs["nrows"] = 0

        header = read_csv(path, **read_kwargs).columns.tolist()

        return self.parse_header(header)

    # ==================================================================================
    # Internal helper methods
    # ==================================================================================
    def parse_header(self, header: Iterable[str]) -> SchemaType:
        """Parse a header and return a schema configured according to it.

        Header fields may contain ``!field!`` or ``?field?`` markers for
        required and optional fields, ``[unit]`` for units, ``(description)``
        for descriptions, and ```validator``` for validator expressions.
        Unmarked fields retain their default ``required`` value.

        Returns
        -------
        SchemaType
            A schema configured according to ``header``.
        """
        schema = self._object.copy(deep=True)
        model = schema._model

        # ------------------------------------------------------------------------------
        if "field" not in model.columns:
            raise ValueError("Schema model must contain a 'field' column")

        # ------------------------------------------------------------------------------
        parsed_header_fields = self._parse_header_fields(header)

        for field, required, unit, validator, description in parsed_header_fields:
            if field not in model["field"].values:
                raise ValueError(f"Field '{field}' not found in schema model")

            if required is not None:
                model.at[field, "required"] = required

            if unit is not None:
                model.at[field, "unit"] = unit

            if validator is not None:
                model.at[field, "validator"] = validator

            if description is not None:
                model.at[field, "description"] = description

        return schema

    # ==================================================================================
    def _parse_header_fields(
        self, header: Iterable[str]
    ) -> Iterator[tuple[str, bool | None, str | None, str | None, str | None]]:
        """Parse a header and yield field information.

        Parameters
        ----------
        header : Iterable[str]
            The header to parse.

        Yields
        ------
        tuple[str, bool | None, str | None, str | None, str | None]
            A tuple of the field name, required, unit, validator, and description.
        """
        for header_field in header:
            field = ""  # TODO: Retreive the field name from the header_field, e.g. "!field! [unit] `validator` (description)" -> "field"

            required = None  # TODO: Retreive the field name from the header_field, e.g. "!field! [unit] `validator` (description)" -> True
            required = None  # TODO: Retreive the field name from the header_field, e.g. "?field? [unit] `validator` (description)" -> False
            required = None  # TODO: Retreive the field name from the header_field, e.g. "field [unit] `validator` (description)" -> None

            unit = None  # TODO: Retreive the field name from the header_field, e.g. "!field! [unit] `validator` (description)" -> "unit"
            unit = None  # TODO: Retreive the field name from the header_field, e.g. "!field! `validator` (description)" -> None

            validator = None  # TODO: Retreive the field name from the header_field, e.g. "!field! `validator` (description)" -> "validator"
            validator = None  # TODO: Retreive the field name from the header_field, e.g. "!field! `validator` (description)" -> None

            description = None  # TODO: Retreive the field name from the header_field, e.g. "!field! `validator` (description)" -> "description"
            description = None  # TODO: Retreive the field name from the header_field, e.g. "!field! `validator`" -> None

            yield field, required, unit, validator, description
