# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING, Any, Self, cast

from pysurveyor.abc.pysurveyor_abc import ABSTRACT_CLASS_ATTRIBUTE
from pysurveyor.backend.backend import DataFrame
from pysurveyor.io.abc.table import Table
from pysurveyor.io.accessors.parser_accessor import ParserAccessor
from pysurveyor.io.loader.schema_parser.schema_parser import SchemaParser
from pysurveyor.schema.utils.validation_model import ValidationMode
from pysurveyor.typing.vars import LoaderType

if TYPE_CHECKING:
    from pysurveyor.data.abc.datatable import DataTable


# ======================================================================================
# Abstract base class for schema objects
# ======================================================================================
class Schema(Table[LoaderType]):
    """Base class for schemas that describe and validate table data.

    A schema stores a tabular model and controls how validation is performed.
    Concrete schema classes provide the appropriate loader, parser, and
    :meth:`validate` implementation for their domain.

    Parameters
    ----------
    model : DataFrame
        Tabular schema model. The model must contain a ``field`` column when
        it is used by :class:`SchemaParser`.
    validation_mode : ValidationMode or str, default=ValidationMode.LAZY
        Validation behavior used by the concrete schema implementation.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    parse = cast(ParserAccessor[SchemaParser], ABSTRACT_CLASS_ATTRIBUTE)

    # ==================================================================================
    def __init__(
        self,
        model: DataFrame,
        validation_mode: ValidationMode | str = ValidationMode.LAZY,
    ) -> None:
        """Initialize a schema from a tabular model."""
        super().__init__(model)
        self.validation_mode = validation_mode

    # ==================================================================================
    @property
    def _model(self) -> DataFrame:
        """Return the schema model used by parsers and validators."""
        return self._data

    # ==================================================================================
    @property
    def validation_mode(self) -> ValidationMode:
        """Return the configured validation mode."""
        return self._validation_mode

    @validation_mode.setter
    def validation_mode(self, new_validation_mode: ValidationMode | str) -> None:
        """Set the validation mode, coercing strings to ``ValidationMode``."""
        if isinstance(new_validation_mode, ValidationMode):
            self._validation_mode = new_validation_mode
        else:
            self._validation_mode = ValidationMode(new_validation_mode)

    # ==================================================================================
    def copy(self, deep: bool = True) -> Self:
        """Return a copy of the schema and its model.

        Parameters
        ----------
        deep : bool, default=True
            Whether to use a deep copy of the underlying model.

        Returns
        -------
        Self
            A schema of the same concrete type.
        """
        return self.__class__(
            self._data.copy(deep=deep),
            **self._copy_kwargs(),
        )

    # ==================================================================================
    def _copy_kwargs(self) -> dict[str, Any]:
        """Return common arguments for a concrete schema copy constructor."""
        return {"validation_mode": self.validation_mode}

    # ==================================================================================
    @abstractmethod
    def validate(self, data: DataTable) -> DataTable:
        """Validate table data against this schema."""
        pass
