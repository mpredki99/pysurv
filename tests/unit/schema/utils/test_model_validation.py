# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Unit tests for schema model validation."""

import pandas as pd
import pytest

from pysurveyor.basic.units import Unit
from pysurveyor.schema.model import models
from pysurveyor.schema.utils.model_validation import validate_model
from pysurveyor.validators.abc.validator import Validator

# ======================================================================================
# Cases and constants for testing model validation
# ======================================================================================
EXPECTED_COLUMNS = {
    "required": bool,
    "unit": Unit,
    "validator": Validator,
    "description": str,
}

MODEL_CASES = {
    "control_points": models.control_points,
    "instruments": models.instruments,
    "measurements": models.measurements,
    "points": models.points,
    "stations": models.stations,
    "targets": models.targets,
}


# ======================================================================================
# Tests for model validation
# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
def test_valid_model(model_name, model):
    """Validate that the valid model remains unchanged after validation."""

    validated = validate_model(
        model,
        EXPECTED_COLUMNS,
        model.index,
        model_name,
    )

    # ----------------------------------------------------------------------------------

    assert validated.equals(model)


# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
def test_standardize_model(model_name, model):
    """Validate that the model is standardized correctly after validation."""

    validated = validate_model(
        model.reset_index(),
        EXPECTED_COLUMNS,
        model.index,
        model_name,
    )

    # ----------------------------------------------------------------------------------

    assert validated.equals(model)


# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
@pytest.mark.parametrize("column", EXPECTED_COLUMNS)
def test_missing_column(model_name, model, column):
    """Validate that a missing column raises a ValueError."""

    invalid_model = model.drop(columns=column)

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="Invalid columns"):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)


# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
def test_missing_index_column(model_name, model):
    """Validate that a missing index column raises a ValueError."""

    invalid_model = model.reset_index().drop(columns="field")

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="Invalid index"):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)


# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
def test_additional_column(model_name, model):
    """Validate that an additional column raises a ValueError."""

    invalid_model = model.assign(extra=None)

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="Invalid columns"):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)


# ======================================================================================
@pytest.mark.parametrize(
    ("model_name", "model", "row"),
    [
        (model_name, model, row)
        for model_name, model in MODEL_CASES.items()
        for row in model.index
    ],
)
def test_missing_row(model_name, model, row):
    """Validate that a missing row raises a ValueError."""

    invalid_model = model.drop(index=row)

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="Invalid index"):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)


# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
def test_additional_row_is_rejected(model_name, model):
    """Validate that an additional row raises a ValueError."""

    extra_row = model.iloc[[0]].copy()
    extra_row.index = ["extra"]
    invalid_model = pd.concat([model, extra_row])

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="Invalid index"):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)


# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
def test_duplicate_index_is_rejected(model_name, model):
    """Validate that a duplicate index raises a ValueError."""

    invalid_model = model.copy()
    invalid_model.index = list(model.index[:-1]) + [model.index[0]]

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="Invalid index"):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)


# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
def test_invalid_index_name_is_rejected(model_name, model):
    invalid_model = model.copy()
    invalid_model.index.name = "invalid"

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="Invalid index"):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)


# ======================================================================================
@pytest.mark.parametrize(("model_name", "model"), MODEL_CASES.items())
@pytest.mark.parametrize(
    ("column", "value"),
    [
        ("required", "not-a-boolean"),
        ("unit", "not-a-unit"),
        ("validator", "not-a-validator"),
        ("description", 123),
    ],
)
def test_invalid_value_type_is_rejected(model_name, model, column, value):
    """Validate that an invalid value type raises a ValueError."""

    invalid_model = model.copy()
    invalid_model[column] = invalid_model[column].astype(object)
    invalid_model[column] = value

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match=column):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)


# ======================================================================================
@pytest.mark.parametrize(
    ("model_name", "model", "row"),
    [
        (model_name, model, row)
        for model_name, model in MODEL_CASES.items()
        for row in model.index
    ],
)
@pytest.mark.parametrize("column", ["unit", "validator", "description"])
def test_empty_value_is_allowed(model_name, model, row, column):
    """Validate that empty metadata values are allowed for any row."""

    valid_model = model.copy()
    valid_model[column] = valid_model[column].astype(object)
    valid_model.at[row, column] = pd.NA

    validated = validate_model(valid_model, EXPECTED_COLUMNS, model.index, model_name)

    # ----------------------------------------------------------------------------------

    assert pd.isna(validated.at[row, column])


# ======================================================================================
@pytest.mark.parametrize(
    ("model_name", "model", "row"),
    [
        (model_name, model, row)
        for model_name, model in MODEL_CASES.items()
        for row in model.index
    ],
)
def test_empty_required(model_name, model, row):
    """Validate that the 'required' column cannot contain an empty value."""

    invalid_model = model.copy()
    invalid_model["required"] = invalid_model["required"].astype(object)
    invalid_model.at[row, "required"] = pd.NA

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="required"):
        validate_model(invalid_model, EXPECTED_COLUMNS, model.index, model_name)
