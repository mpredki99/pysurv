# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from typing import Any, cast

from pysurveyor.abc.pysurveyor_abc import ABSTRACT_CLASS_ATTRIBUTE
from pysurveyor.backend.backend import DataFrame, Index
from pysurveyor.basic.units import Unit
from pysurveyor.schema.abc.schema import Schema
from pysurveyor.schema.utils.model_validation import (
    assert_valid_columns_dtype,
    assert_valid_columns_existence,
    assert_valid_index,
    validate_model,
)
from pysurveyor.schema.utils.validation_mode import ValidationMode
from pysurveyor.typing.vars import LoaderType
from pysurveyor.validators.abc.validator import Validator


# ======================================================================================
# Schema with a fixed model structure and declared column types.
# ======================================================================================
class StrictSchema(Schema[LoaderType]):
    """``Schema`` with a fixed model structure and declared column types.

    ``StrictSchema`` enforces the declared model structure and validates non-null column
    values against their expected types. Model assignments are validated before being
    committed, ensuring that failed assignments leave the existing model unchanged.

    Subclasses must define the expected index and may override the expected column
    definitions to specify their schema requirements. Once initialized, schema
    attributes cannot be modified or deleted directly.
    """

    __pysurveyor_abstract__ = True

    # ----------------------------------------------------------------------------------
    _expected_columns: dict[str, type] = {
        "required": bool,
        "unit": Unit,
        "validator": Validator,
        "description": str,
    }
    # ----------------------------------------------------------------------------------
    _expected_index = cast(Index, ABSTRACT_CLASS_ATTRIBUTE)

    # ==================================================================================
    # Dunder methods
    # ==================================================================================
    def __init__(
        self,
        model: DataFrame,
        validation_mode: ValidationMode | str = ValidationMode.LAZY,
    ) -> None:

        super().__init__(model, validation_mode=validation_mode)
        # ------------------------------------------------------------------------------
        self._assert_valid_model()
        # ------------------------------------------------------------------------------
        self._is_initialized = True

    # ==================================================================================
    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, "_is_initialized", False):
            raise AttributeError(
                f"'{type(self).__name__}' does not allow modifying attribute {name!r}."
            )

        object.__setattr__(self, name, value)

    # ==================================================================================
    # Attribute protection
    # ==================================================================================
    def __delattr__(self, name: str) -> None:
        raise AttributeError(
            f"'{type(self).__name__}' does not allow deleting attribute {name!r}."
        )

    # ==================================================================================
    def __getitem__(self, key: Any) -> Any:
        return self._model[key]

    # ==================================================================================
    # Data assignment validation
    # ==================================================================================
    def __setitem__(self, key: Any, value: Any) -> None:
        """Apply assignment, then validate and commit."""
        schema_name = type(self).__name__

        # Copy the model to avoid modifying the original until validation passes -------
        candidate = self._model.copy(deep=True)
        candidate[key] = value

        # ------------------------------------------------------------------------------
        assert_valid_index(
            candidate.index,
            self._expected_index,
            schema_name,
        )
        # ------------------------------------------------------------------------------
        assert_valid_columns_existence(
            candidate.columns,
            self._expected_columns.keys(),
            schema_name,
        )
        # ------------------------------------------------------------------------------
        assert_valid_columns_dtype(
            candidate,
            self._expected_columns,
            schema_name,
        )

        # Commit only after every validation succeeds ----------------------------------
        object.__setattr__(self, "_model", candidate)

    # ==================================================================================
    # Internal methods
    # ==================================================================================
    def _assert_valid_model(self) -> None:
        self._model = validate_model(
            self._model,
            self._expected_columns,
            self._expected_index,
            self.__class__.__name__,
        )
