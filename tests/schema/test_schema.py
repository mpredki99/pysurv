# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pandas as pd
import pytest

from pysurv.exceptions._exceptions import ValidationError
from pysurv.schema import FlexibleSchema, StrictSchema
from pysurv.schema.utils import SchemaValidationMode


# --------------------------------------------------------------------------------------
# Test validation mode setter
# --------------------------------------------------------------------------------------
def test_schema_validation_mode():
    schema = StrictSchema(pd.DataFrame(), validation_mode=SchemaValidationMode.DISABLED)
    assert schema.validation_mode is SchemaValidationMode.DISABLED

    schema.validation_mode = SchemaValidationMode.EAGER
    assert schema.validation_mode is SchemaValidationMode.EAGER


def test_schema_validation_from_string():
    schema = StrictSchema(pd.DataFrame(), validation_mode="fields")
    assert schema.validation_mode is SchemaValidationMode.FIELDS

    schema.validation_mode = "coerce"
    assert schema.validation_mode is SchemaValidationMode.COERCE


# --------------------------------------------------------------------------------------
# Test validation on DISABLED mode
# --------------------------------------------------------------------------------------
def test_validation_disabled_valid_data(valid_model_file_path, valid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="disabled")

    validated = schema.validate(valid_data.copy())
    assert validated.equals(valid_data)


def test_validation_disabled_invalid_data(valid_model_file_path, invalid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="disabled")

    validated = schema.validate(invalid_data.copy())
    assert validated.equals(invalid_data)


def test_validation_disabled_missing_column_data(
    valid_model_file_path, missing_column_data
):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="disabled")

    validated = schema.validate(missing_column_data.copy())
    assert validated.equals(missing_column_data)


# --------------------------------------------------------------------------------------
# Test validation on FIELDS mode
# --------------------------------------------------------------------------------------
def test_validation_fields_valid_data(valid_model_file_path, valid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="fields")

    validated = schema.validate(valid_data.copy())
    assert validated.equals(valid_data)


def test_validation_fields_invalid_data(valid_model_file_path, invalid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="fields")

    validated = schema.validate(invalid_data.copy())
    assert validated.equals(invalid_data)


def test_validation_fields_missing_column_data(
    valid_model_file_path, missing_column_data
):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="fields")

    with pytest.raises(ValidationError):
        schema.validate(missing_column_data.copy())


# --------------------------------------------------------------------------------------
# Test validation on EAGER mode
# --------------------------------------------------------------------------------------
def test_validation_eager_valid_data(valid_model_file_path, valid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="eager")

    validated = schema.validate(valid_data.copy())
    assert validated["unvalidated"].equals(valid_data["unvalidated"])
    # Assert field were converted to correct type
    assert validated["text"].dtype == "str"
    assert validated["distance"].dtype == "float64"
    assert validated["angle"].dtype == "float64"


def test_validation_eager_invalid_data(valid_model_file_path, invalid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="eager")

    with pytest.raises(ValidationError):
        schema.validate(invalid_data.copy())


def test_validation_eager_missing_column_data(
    valid_model_file_path, missing_column_data
):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="eager")

    with pytest.raises(ValidationError):
        schema.validate(missing_column_data.copy())


# --------------------------------------------------------------------------------------
# Test validation on LAZY mode
# --------------------------------------------------------------------------------------
def test_validation_lazy_valid_data(valid_model_file_path, valid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="lazy")

    validated = schema.validate(valid_data.copy())
    assert validated["unvalidated"].equals(valid_data["unvalidated"])
    # Assert field were converted to correct type
    assert validated["text"].dtype == "str"
    assert validated["distance"].dtype == "float64"
    assert validated["angle"].dtype == "float64"


def test_validation_lazy_invalid_data(valid_model_file_path, invalid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="lazy")

    with pytest.raises(ValidationError):
        schema.validate(invalid_data.copy())


def test_validation_lazy_missing_column_data(
    valid_model_file_path, missing_column_data
):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="lazy")

    with pytest.raises(ValidationError):
        schema.validate(missing_column_data.copy())


# --------------------------------------------------------------------------------------
# Test validation on COERCE mode
# --------------------------------------------------------------------------------------
def test_validation_coerce_valid_data(valid_model_file_path, valid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="coerce")

    validated = schema.validate(valid_data.copy())
    assert validated["unvalidated"].equals(valid_data["unvalidated"])
    # Assert field were converted to correct type
    assert validated["text"].dtype == "str"
    assert validated["distance"].dtype == "float64"
    assert validated["angle"].dtype == "float64"


def test_validation_coerce_invalid_data(valid_model_file_path, invalid_data):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="coerce")

    validated = schema.validate(invalid_data.copy())
    assert validated["unvalidated"].equals(invalid_data["unvalidated"])
    # Assert field were converted to correct type
    assert validated["text"].dtype == "str"
    assert validated["distance"].dtype == "float64"
    assert validated["angle"].dtype == "float64"
    # Assert invalid values were replaced with empty values
    assert validated.isna().sum().sum() > invalid_data.isna().sum().sum()


def test_validation_coerce_missing_column_data(
    valid_model_file_path, missing_column_data
):
    schema = FlexibleSchema.from_csv(valid_model_file_path, validation_mode="coerce")

    with pytest.raises(ValidationError):
        schema.validate(missing_column_data.copy())
