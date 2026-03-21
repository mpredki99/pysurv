# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.


import pytest

from pysurv.schema import MeasurementsSchema, StrictSchema
from pysurv.typing import AngleUnit, DistanceUnit, PySurvUnit
from pysurv.validators import Between, PySurvValidator


def test_initial_type(invalid_model_file_path):
    with pytest.raises(ValueError):
        StrictSchema.from_csv(invalid_model_file_path)


# --------------------------------------------------------------------------------------
# Test setting value on __getitem__ method
# --------------------------------------------------------------------------------------
def test_assign_valid_column():
    schema = MeasurementsSchema()

    schema["required"] = True
    assert (schema["required"] == True).all()

    schema["unit"] = DistanceUnit.METERS
    assert (schema["unit"].apply(lambda x: isinstance(x, PySurvUnit))).all()

    schema["validator"] = Between(0, 100)
    assert (schema["unit"].apply(lambda x: isinstance(x, PySurvUnit))).all()

    schema["description"] = "New description"
    assert (schema["description"].apply(lambda x: isinstance(x, str))).all()


def test_assign_invalid_column():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema["required"] = None

    with pytest.raises(TypeError):
        schema["unit"] = "Invalid type"

    with pytest.raises(TypeError):
        schema["validator"] = "Invalid type"

    with pytest.raises(TypeError):
        schema["description"] = None


def test_add_new_column():
    schema = MeasurementsSchema()

    with pytest.raises(KeyError):
        schema["new"] = None


# --------------------------------------------------------------------------------------
# Test setting value on Loc indexer
# --------------------------------------------------------------------------------------
def test_assign_valid_value_loc():
    schema = MeasurementsSchema()

    schema.loc["sd", "required"] = True
    assert schema.loc["sd", "required"] == True

    schema.loc["sd", "unit"] = DistanceUnit.METERS
    assert isinstance(schema.loc["sd", "unit"], PySurvUnit)

    schema.loc["sd", "validator"] = Between(0, 100)
    assert isinstance(schema.loc["sd", "validator"], PySurvValidator)

    schema.loc["sd", "description"] = "New description"
    assert isinstance(schema.loc["sd", "description"], str)


def test_assign_valid_column_loc():
    schema = MeasurementsSchema()

    schema.loc[:, "required"] = True
    assert (schema.loc[:, "required"] == True).all()

    schema.loc[:, "unit"] = DistanceUnit.METERS
    assert (schema.loc[:, "unit"].apply(lambda x: isinstance(x, PySurvUnit))).all()

    schema.loc[:, "validator"] = Between(0, 100)
    assert (
        schema.loc[:, "validator"].apply(lambda x: isinstance(x, PySurvValidator))
    ).all()

    schema.loc[:, "description"] = "New description"
    assert (schema.loc[:, "description"].apply(lambda x: isinstance(x, str))).all()


def test_assign_valid_row_loc():
    schema = MeasurementsSchema()

    schema.loc["sd"] = True, AngleUnit.GRAD, Between(0, 100), "New description"
    assert schema.loc["sd", "required"] == True
    assert isinstance(schema.loc["sd", "unit"], PySurvUnit)
    assert isinstance(schema.loc["sd", "validator"], PySurvValidator)
    assert isinstance(schema.loc["sd", "description"], str)


def test_assign_invalid_value_loc():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema.loc["sd", "required"] = None

    with pytest.raises(TypeError):
        schema.loc["sd", "unit"] = "Invalid type"

    with pytest.raises(TypeError):
        schema.loc["sd", "validator"] = "Invalid type"

    with pytest.raises(TypeError):
        schema.loc["sd", "description"] = None


def test_assign_invalid_column_loc():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema.loc[:, "required"] = None

    with pytest.raises(TypeError):
        schema.loc[:, "unit"] = "Invalid type"

    with pytest.raises(TypeError):
        schema.loc[:, "validator"] = "Invalid type"

    with pytest.raises(TypeError):
        schema.loc[:, "description"] = None


def test_assign_invalid_row_loc():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema.loc["sd"] = None, AngleUnit.GRAD, "Invalid type", "New description"


def test_add_new_column_loc():
    schema = MeasurementsSchema()

    with pytest.raises(KeyError):
        schema.loc[:, "new"] = None


def test_add_new_row_loc():
    schema = MeasurementsSchema()

    with pytest.raises(KeyError):
        schema.loc["new", :] = None


def test_add_new_item_loc():
    schema = MeasurementsSchema()

    with pytest.raises(KeyError):
        schema.loc["new_row", "new_col"] = None


# --------------------------------------------------------------------------------------
# Test setting value on Iloc indexer
# --------------------------------------------------------------------------------------
def test_assign_valid_value_iloc():
    schema = MeasurementsSchema()

    schema.iloc[10, 0] = True
    assert schema.iloc[10, 0] == True

    schema.iloc[10, 1] = DistanceUnit.METERS
    assert isinstance(schema.iloc[10, 1], PySurvUnit)

    schema.iloc[10, 2] = Between(0, 100)
    assert isinstance(schema.iloc[10, 2], PySurvValidator)

    schema.iloc[10, 3] = "New description"
    assert isinstance(schema.iloc[10, 3], str)


def test_assign_valid_column_iloc():
    schema = MeasurementsSchema()

    schema.iloc[:, 0] = True
    assert (schema.iloc[:, 0] == True).all()

    schema.iloc[:, 1] = DistanceUnit.METERS
    assert (schema.iloc[:, 1].apply(lambda x: isinstance(x, PySurvUnit))).all()

    schema.iloc[:, 2] = Between(0, 100)
    assert (schema.iloc[:, 2].apply(lambda x: isinstance(x, PySurvValidator))).all()

    schema.iloc[:, 3] = "New description"
    assert (schema.iloc[:, 3].apply(lambda x: isinstance(x, str))).all()


def test_assign_valid_row_iloc():
    schema = MeasurementsSchema()

    schema.iloc[10] = True, AngleUnit.GRAD, Between(0, 100), "New description"
    assert schema.iloc[10, 0] == True
    assert isinstance(schema.iloc[10, 1], PySurvUnit)
    assert isinstance(schema.iloc[10, 2], PySurvValidator)
    assert isinstance(schema.iloc[10, 3], str)


def test_assign_invalid_value_iloc():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema.iloc[10, 0] = None

    with pytest.raises(TypeError):
        schema.iloc[10, 1] = "Invalid type"

    with pytest.raises(TypeError):
        schema.iloc[10, 2] = "Invalid type"

    with pytest.raises(TypeError):
        schema.iloc[10, 3] = None


def test_assign_invalid_column_iloc():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema.iloc[:, 0] = None

    with pytest.raises(TypeError):
        schema.iloc[:, 1] = "Invalid type"

    with pytest.raises(TypeError):
        schema.iloc[:, 2] = "Invalid type"

    with pytest.raises(TypeError):
        schema.iloc[:, 3] = None


def test_assign_invalid_row_iloc():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema.iloc[10] = None, AngleUnit.GRAD, "Invalid type", "New description"


def test_add_new_column_iloc():
    schema = MeasurementsSchema()

    with pytest.raises(IndexError):
        schema.iloc[:, 4] = None


def test_add_new_row_iloc():
    schema = MeasurementsSchema()

    with pytest.raises(IndexError):
        schema.iloc[27, :] = None


def test_add_new_item_iloc():
    schema = MeasurementsSchema()

    with pytest.raises(IndexError):
        schema.iloc[27, 4] = None


# --------------------------------------------------------------------------------------
# Test setting value on At indexer
# --------------------------------------------------------------------------------------
def test_assign_valid_value_at():
    schema = MeasurementsSchema()

    schema.at["sd", "required"] = True
    assert schema.at["sd", "required"] == True

    schema.at["sd", "unit"] = DistanceUnit.METERS
    assert isinstance(schema.at["sd", "unit"], PySurvUnit)

    schema.at["sd", "validator"] = Between(0, 100)
    assert isinstance(schema.at["sd", "validator"], PySurvValidator)

    schema.at["sd", "description"] = "New description"
    assert isinstance(schema.at["sd", "description"], str)


def test_assign_invalid_value_at():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema.at["sd", "required"] = None

    with pytest.raises(TypeError):
        schema.at["sd", "unit"] = "Invalid type"

    with pytest.raises(TypeError):
        schema.at["sd", "validator"] = "Invalid type"

    with pytest.raises(TypeError):
        schema.at["sd", "description"] = None


def test_add_new_item_at():
    schema = MeasurementsSchema()

    with pytest.raises(KeyError):
        schema.at["new_row", "new_col"] = None


# --------------------------------------------------------------------------------------
# Test setting value on Iat indexer
# --------------------------------------------------------------------------------------
def test_assign_valid_value_iat():
    schema = MeasurementsSchema()

    schema.iat[10, 0] = True
    assert schema.iat[10, 0] == True

    schema.iat[10, 1] = DistanceUnit.METERS
    assert isinstance(schema.iat[10, 1], PySurvUnit)

    schema.iat[10, 2] = Between(0, 100)
    assert isinstance(schema.iat[10, 2], PySurvValidator)

    schema.iat[10, 3] = "New description"
    assert isinstance(schema.iat[10, 3], str)


def test_assign_invalid_value_iat():
    schema = MeasurementsSchema()

    with pytest.raises(TypeError):
        schema.iat[10, 0] = None

    with pytest.raises(TypeError):
        schema.iat[10, 1] = "Invalid type"

    with pytest.raises(TypeError):
        schema.iat[10, 2] = "Invalid type"

    with pytest.raises(TypeError):
        schema.iat[10, 3] = None


def test_add_new_item_iat():
    schema = MeasurementsSchema()

    with pytest.raises(IndexError):
        schema.iat[27, 4] = None
