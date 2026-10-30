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
from pysurveyor.basic.units import Unit, parse_unit
from pysurveyor.io.config.abc.config import IOConfig
from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.loader.abc.object_bound_loader import ObjectBoundLoader
from pysurveyor.io.source.abc.source import Source
from pysurveyor.io.source.file_sources import CSVSource
from pysurveyor.io.source.resolver import resolve_source
from pysurveyor.typing.vars import SchemaType
from pysurveyor.validators.abc.validator import Validator
from pysurveyor.validators.utils.parser import parse_validator

# ======================================================================================
# Header patterns compiled once
# ======================================================================================
_HEADER_PATTERNS = {
    "required": re.compile(r"!([^!]+)!"),
    "optional": re.compile(r"\?([^?]+)\?"),
    "unit": re.compile(r"\[([^\]]+)\]"),
    "validator": re.compile(r"`([^`]*)`"),
    "description": re.compile(r"\(([^()]*)\)"),
}


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
        if model.index.name != "field":
            raise ValueError("Schema model must have index named 'field'")

        # ------------------------------------------------------------------------------
        parsed_header_fields = self._parse_header_fields(header)

        for field, required, unit, validator, description in parsed_header_fields:

            if field not in model.index:
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
    ) -> Iterator[tuple[str, bool | None, Unit | None, Validator | None, str | None]]:
        """Parse header fields.

        Supported annotations
        ---------------------
        !field!       Required field
        ?field?       Optional field
        [unit]        Unit
        `validator`   Validator expression
        (description) Field description

        Annotations may appear in any order, with arbitrary whitespace.
        """
        for header_field in header:
            remaining = header_field.strip()

            # --------------------------------------------------------------------------
            # Extract the required/optional marker
            # --------------------------------------------------------------------------
            remaining, required_match = self._extract_header_component(
                remaining,
                _HEADER_PATTERNS["required"],
                "required marker",
                preserve_value=True,
            )
            remaining, optional_match = self._extract_header_component(
                remaining,
                _HEADER_PATTERNS["optional"],
                "optional marker",
                preserve_value=True,
            )

            if required_match is not None and optional_match is not None:
                raise ValueError(
                    "Header field cannot be both required and optional: "
                    f"{header_field!r}"
                )

            required = (
                True
                if required_match is not None
                else False if optional_match is not None else None
            )

            # --------------------------------------------------------------------------
            # Extract the unit
            # --------------------------------------------------------------------------
            remaining, unit = self._extract_header_component(
                remaining, _HEADER_PATTERNS["unit"], "unit"
            )
            unit = parse_unit(unit.strip()) if unit is not None else None

            # --------------------------------------------------------------------------
            # Extract the validator
            # --------------------------------------------------------------------------
            remaining, validator = self._extract_header_component(
                remaining, _HEADER_PATTERNS["validator"], "validator"
            )
            validator = (
                parse_validator(validator.strip()) if validator is not None else None
            )

            # --------------------------------------------------------------------------
            # Extract the description
            # --------------------------------------------------------------------------
            remaining, description = self._extract_header_component(
                remaining, _HEADER_PATTERNS["description"], "description"
            )

            # --------------------------------------------------------------------------
            # Remove whitespace and separating annotations to get the field name
            # --------------------------------------------------------------------------
            field = remaining.strip().strip(",;").strip()

            if not field:
                raise ValueError(f"Missing field name in header: {header_field!r}")

            # --------------------------------------------------------------------------
            yield field, required, unit, validator, description

    # ==================================================================================
    @staticmethod
    def _extract_header_component(
        text: str,
        pattern: re.Pattern[str],
        component: str,
        preserve_value: bool = False,
    ) -> tuple[str, str | None]:
        """Extract a single annotation and return the remaining text.

        If preserve_value is True, retain the captured value in the text.

        Raises
        ------
        ValueError
            If multiple annotations of the same kind occur.
        """
        match = pattern.search(text)

        # ------------------------------------------------------------------------------
        if match is None:
            return text, None

        # ------------------------------------------------------------------------------
        if pattern.search(text, match.end()) is not None:
            raise ValueError(f"Multiple {component} annotations in header: {text!r}")

        # ------------------------------------------------------------------------------
        value = match.group(1).strip()
        replacement = value if preserve_value else ""
        remaining = text[: match.start()] + replacement + text[match.end() :]

        return remaining, value
