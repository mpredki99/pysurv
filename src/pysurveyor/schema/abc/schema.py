# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING, Any, Self, cast

from pysurveyor.abc.pysurveyor_abc import ABSTRACT_CLASS_ATTRIBUTE
from pysurveyor.backend.backend import NA, DataFrame, Index, isna
from pysurveyor.exceptions.exceptions import ValidationError
from pysurveyor.io.abc.table import Table
from pysurveyor.io.accessors.parser_accessor import ParserAccessor
from pysurveyor.io.loader.schema_parser.schema_parser import SchemaParser
from pysurveyor.schema.utils.validation_mode import ValidationMode
from pysurveyor.typing.vars import LoaderType
from pysurveyor.validators.utils.validator_mode import ValidatorMode

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
    parse = cast(ParserAccessor[SchemaParser[Self]], ABSTRACT_CLASS_ATTRIBUTE)

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
    # Properties
    # ==================================================================================
    @property
    def _model(self) -> DataFrame:
        """Return the schema model."""
        return self._data

    @_model.setter
    def _model(self, value: DataFrame) -> None:
        """Set the schema model."""
        self._data = value

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
    @property
    def data_columns(self) -> Index:
        """Return the data column names."""
        return self._model.index

    # ==================================================================================
    # Copy methods
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
    # Public validation method
    # ==================================================================================
    def validate(self, data: DataTable) -> DataTable:
        """Validate table data against this schema and validation mode."""
        mode_methods = {
            ValidationMode.FIELDS: self._validate_fields,
            ValidationMode.EAGER: self._validate_eager,
            ValidationMode.LAZY: self._validate_lazy,
            ValidationMode.EMPTY: self._validate_empty,
            ValidationMode.DISABLED: None,
        }
        validation_method = mode_methods[self.validation_mode]

        # ------------------------------------------------------------------------------
        if validation_method is None:
            return data

        # ------------------------------------------------------------------------------
        return validation_method(data)

    # ==================================================================================
    # Internal validation methods for each mode
    # ==================================================================================
    def _validate_fields(self, data: DataTable) -> DataTable:
        """Check only existence of required fields."""
        required_mask = self._model["required"].fillna(False).astype(bool)
        required_fields = set(self._model.index[required_mask])
        data_fields = set(data.columns)

        # ------------------------------------------------------------------------------
        if required_fields.issubset(data_fields):
            return data

        # ------------------------------------------------------------------------------
        raise ValidationError(
            f"Missing mandatory columns: {required_fields - data_fields}"
        )

    # ==================================================================================
    def _validate_eager(self, data: DataTable) -> DataTable:
        """Validate data. Raise validation error on first invalid column."""
        data = self._validate_fields(data)

        # ------------------------------------------------------------------------------
        validators = self._model["validator"]

        # ------------------------------------------------------------------------------
        for field in data.columns:
            if field not in validators.index:
                continue

            validator = validators[field]

            # --------------------------------------------------------------------------
            if isna(validator):
                continue

            # Ensure validator is in RAISE mode for eager validation -------------------
            validator = deepcopy(validator)
            validator.mode = ValidatorMode.RAISE

            # --------------------------------------------------------------------------
            _, data[field] = validator(data[field])
        # ------------------------------------------------------------------------------
        return data

    # ==================================================================================
    def _validate_lazy(self, data: DataTable) -> DataTable:
        """Validate data. Raise validation error with summary of all invalid data."""
        data = self._validate_fields(data)

        # ------------------------------------------------------------------------------
        invalid, data = self._evaluate_data(data)

        # ------------------------------------------------------------------------------
        if not invalid.empty:
            # Show only rows that have at least one invalid value ----------------------
            invalid_summary = invalid[invalid.any(axis=1)]
            raise ValidationError("Invalid values found:\n" f"{invalid_summary}")

        # ------------------------------------------------------------------------------
        return data

    # ==================================================================================
    def _validate_empty(self, data: DataTable) -> DataTable:
        """Validate data. Replace invalid data with empty values."""
        data = self._validate_fields(data)

        # ------------------------------------------------------------------------------
        invalid, data = self._evaluate_data(data)

        # ------------------------------------------------------------------------------
        if invalid.any().any():
            data[invalid] = NA

        # ------------------------------------------------------------------------------
        return data

    # ==================================================================================
    # Validation helper method
    # ==================================================================================
    def _evaluate_data(self, data: DataTable) -> tuple[DataFrame, DataTable]:
        """Validate data. Accumulate invalid values and return it with the data."""
        invalid_masks = {}
        validators = self._model["validator"]

        # ------------------------------------------------------------------------------
        for field in data.columns:
            if field not in validators.index:
                continue

            validator = validators[field]

            # --------------------------------------------------------------------------
            if isna(validator):
                continue

            # Avoid mutating validator object stored in schema model -------------------
            validator = deepcopy(validator)
            validator.mode = ValidatorMode.RETURN

            # --------------------------------------------------------------------------
            mask, data[field] = validator(data[field])

            # --------------------------------------------------------------------------
            if not mask.all():
                invalid_masks[field] = ~mask

        # ------------------------------------------------------------------------------
        if not invalid_masks:
            return DataFrame(index=data.index), data

        invalid = DataFrame(invalid_masks, index=data.index)
        return invalid, data
