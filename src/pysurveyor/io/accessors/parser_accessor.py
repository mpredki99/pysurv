# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Generic, cast

from pysurveyor.abc.pysurveyor_abc import is_pysurveyor_abstract
from pysurveyor.typing.vars import SchemaParserType

if TYPE_CHECKING:
    from pysurveyor.schema.abc.schema import Schema


# ======================================================================================
# Class descriptor for object-bound schema parsers
# ======================================================================================
class ParserAccessor(Generic[SchemaParserType]):
    """Expose a schema parser bound to a schema object.

    Instance access caches the parser on the schema instance. Class access
    creates a parser bound to a fresh default schema instance without
    caching the parser on the class.

    Attributes
    ----------
    _constructor : type[SchemaParserType]
        Parser class used to create the bound parser instance.
    _attr : str
        Name of the instance attribute used to cache the parser.
    """

    # ==================================================================================
    def __init__(self, constructor: type[SchemaParserType]) -> None:
        """Initialize the accessor with the parser constructor."""
        self._constructor = constructor
        self._attr = ""

    # ==================================================================================
    def __set_name__(self, owner: type[Any], name: str) -> None:
        """Register the instance attribute used to cache the parser."""
        self._attr = f"_schema_parser_of_{name}"

    # ==================================================================================
    def __get__(
        self,
        obj: Schema[Any] | None,
        owner: type[Schema[Any]],
    ) -> SchemaParserType:
        """Return the parser bound to a schema instance.

        Class access creates a parser bound to a new schema instance.
        Instance access creates and caches a parser bound to the accessed
        schema object.

        Parameters
        ----------
        obj : Schema | None
            Schema instance through which the descriptor is accessed.
        owner : type[Schema]
            Class that owns the parser descriptor.

        Returns
        -------
        SchemaParserType
            Parser instance constructed with a schema object.

        Raises
        ------
        TypeError
            If class access is attempted on an abstract schema class or
            the owning schema class cannot be instantiated without arguments.
        """
        if obj is None:
            return self._on_class_access(owner)

        if self._attr not in obj.__dict__:
            setattr(obj, self._attr, self._constructor(obj))

        return cast(SchemaParserType, obj.__dict__[self._attr])

    # ==================================================================================
    def _on_class_access(self, owner: type[Schema[Any]]) -> SchemaParserType:
        """Create a parser bound to a fresh schema instance."""
        if is_pysurveyor_abstract(owner):
            raise TypeError(
                f"Cannot use {owner.__name__}.parse on an abstract schema class. "
                f"Access parse on a concrete schema class or on a schema instance."
            )
        # ------------------------------------------------------------------------------
        try:
            schema = cast(Any, owner)()
        except TypeError as exc:
            raise TypeError(
                f"Cannot use {owner.__name__}.parse on the class because "
                f"{owner.__name__} could not be instantiated to bind the parser. "
                f"Access parse on an existing schema instance instead."
            ) from exc

        return self._constructor(schema)
