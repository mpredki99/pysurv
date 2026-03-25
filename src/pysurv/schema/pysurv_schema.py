# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from abc import ABC
from copy import deepcopy
from typing import Any, Tuple

import pandas as pd

from pysurv.configs.format_config import CSVConfig
from pysurv.data.pysurv_table import PySurvTable
from pysurv.exceptions._exceptions import ValidationError
from pysurv.validators.utils import PySurvValidatorMode

from .field_compilers import compile_unit, compile_validator
from .utils import SchemaValidationMode


class PySurvSchema(ABC):
    """
    Abstract base class for PySurv schema.
    Handles schema model DataFrame creation, validation logic,
    and defines public interface.
    """

    def __init__(
        self,
        model: pd.DataFrame,
        validation_mode: SchemaValidationMode | str = SchemaValidationMode.LAZY,
    ) -> None:
        super().__init__()
        self.validation_mode = validation_mode
        self._model = self._normalize_model(model)

    # ----------------------------------------------------------------------------------
    # Factory methods
    # ----------------------------------------------------------------------------------
    @classmethod
    def from_csv(
        cls,
        path: str,
        config: CSVConfig | None = None,
        validation_mode: SchemaValidationMode | str = SchemaValidationMode.LAZY,
    ) -> "PySurvSchema":
        config = config or CSVConfig()
        model = pd.read_csv(path, **config.kwargs)
        return cls(model, validation_mode=validation_mode)

    # ----------------------------------------------------------------------------------
    # Dunder methods
    # ----------------------------------------------------------------------------------
    def __getattr__(self, name: str) -> Any:
        if name == "_model":
            raise AttributeError("_model not initialized")

        try:
            return getattr(self._model, name)
        except AttributeError:
            raise AttributeError(
                f"{self.__class__.__name__} does not have attribute {name}"
            )

    def __repr__(self) -> str:
        repr_model = self._model.copy()
        repr_model = self._repr_field(repr_model, "unit")
        repr_model = self._repr_field(repr_model, "validator")

        return (
            f"validation_mode: {self.validation_mode!r}"
            "\n"
            f"{repr_model.to_string(max_colwidth=None)}"
        )

    def __str__(self) -> str:
        str_model = self._model.copy()
        str_model = self._str_field(str_model, "validator")

        return f"validation_mode: {self.validation_mode}" "\n" f"{str_model}"

    # ----------------------------------------------------------------------------------
    # Properties
    # ----------------------------------------------------------------------------------
    @property
    def validation_mode(self) -> SchemaValidationMode:
        return self._validation_mode

    @validation_mode.setter
    def validation_mode(self, value: SchemaValidationMode | str) -> None:
        if isinstance(value, SchemaValidationMode):
            self._validation_mode = value
        else:
            self._validation_mode = SchemaValidationMode(value)

    # ----------------------------------------------------------------------------------
    # Public interface
    # ----------------------------------------------------------------------------------
    def to_data_frame(self) -> pd.DataFrame:
        return self._model.copy(deep=True)

    def validate(self, data: PySurvTable) -> PySurvTable:
        """Validate data according to schema and validation mode."""
        modes = {
            SchemaValidationMode.FIELDS: self._validate_fields,
            SchemaValidationMode.EAGER: self._validate_eager,
            SchemaValidationMode.LAZY: self._validate_lazy,
            SchemaValidationMode.COERCE: self._validate_coerce,
            SchemaValidationMode.DISABLED: None,
        }
        validation_method = modes[self.validation_mode]

        if validation_method is None:
            return data

        return validation_method(data)

    # ----------------------------------------------------------------------------------
    # Internal validators
    # ----------------------------------------------------------------------------------
    def _validate_fields(self, data: PySurvTable) -> PySurvTable:
        """Check only existence of required fields."""
        required_fields = set(self._model[self._model["required"] == True].index)
        data_fields = set(data.columns)

        if required_fields.issubset(data.columns):
            return data

        raise ValidationError(
            f"Missing mandatory columns: {required_fields - data_fields}"
        )

    def _validate_eager(self, data: PySurvTable) -> PySurvTable:
        """Validate data. Raise validation error on first invalid column."""
        data = self._validate_fields(data)

        validators = self._model["validator"]

        for field in data.columns:
            validator = validators[field]
            if pd.isna(validator):
                continue

            _, data[field] = validator(data[field])

        return data

    def _evaluate_data(self, data: PySurvTable) -> Tuple[pd.DataFrame, PySurvTable]:
        """Validate data. Accumulate invalid values and return it with the data."""
        invalid = pd.DataFrame(index=data.index)
        validators = self._model["validator"]

        for field in data.columns:
            # Avoid mutating validator object stored in schema model
            validator = deepcopy(validators[field])
            if pd.isna(validator):
                continue

            validator.mode = PySurvValidatorMode.RETURN
            mask, data[field] = validator(data[field])

            if not mask.all():
                invalid[field] = ~mask

        return invalid, data

    def _validate_lazy(self, data: PySurvTable) -> PySurvTable:
        """Validate data. Raise validation error with summary of all invalid data."""
        data = self._validate_fields(data)

        invalid, data = self._evaluate_data(data)

        if invalid.any().any():
            raise ValidationError("Invalid values found" "\n" f"{invalid}")

        return data

    def _validate_coerce(self, data: PySurvTable) -> PySurvTable:
        """Validate data. Replace invalid data with empty values."""
        data = self._validate_fields(data)

        invalid, data = self._evaluate_data(data)

        if invalid.any().any():
            data[invalid] = pd.NA

        return data

    # ----------------------------------------------------------------------------------
    # Internal model transformers
    # ----------------------------------------------------------------------------------
    def _normalize_model(self, model: pd.DataFrame) -> pd.DataFrame:
        """Normalize and parse the raw schema model DataFrame."""
        if "field" in model:
            model = model.set_index("field")
        return self._parse_columns(model)

    def _parse_columns(self, model: pd.DataFrame) -> pd.DataFrame:
        """Parse input data into specific dtype."""
        model = self._parse_units(model) if "unit" in model else model
        model = self._parse_validators(model) if "validator" in model else model
        return model

    def _parse_units(self, model: pd.DataFrame) -> pd.DataFrame:
        """Replace unit strings with PySurvUnit objects."""
        model["unit"] = [compile_unit(value) for value in model["unit"]]
        return model

    def _parse_validators(self, model: pd.DataFrame) -> pd.DataFrame:
        """Replace validator strings with PySurvValidator objects."""
        model["validator"] = [compile_validator(value) for value in model["validator"]]
        return model

    def _repr_field(self, repr_model: pd.DataFrame, field: str) -> pd.DataFrame:
        """Replace field's object with their dev representation."""
        repr_model[field] = [repr(v) if v is not None else v for v in repr_model[field]]
        return repr_model

    def _str_field(self, repr_model: pd.DataFrame, field: str) -> pd.DataFrame:
        """Replace field's object with thier string representation."""
        repr_model[field] = [str(v) if v is not None else v for v in repr_model[field]]
        return repr_model
