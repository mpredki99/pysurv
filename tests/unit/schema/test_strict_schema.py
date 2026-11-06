# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests for strict schema behavior."""

import pandas as pd
import pytest

from pysurveyor.basic import DistanceUnit
from pysurveyor.schema import (
    ControlPointsSchema,
    InstrumentsSchema,
    MeasurementsSchema,
    PointsSchema,
    StationsSchema,
    TargetsSchema,
    ValidationMode,
)
from pysurveyor.validators.conversion_validators import IsNumeric

# ======================================================================================
# Strict schema classes for testing
# ======================================================================================
SCHEMA_CLASSES_AND_INDEXES = {
    ControlPointsSchema: ControlPointsSchema().index,
    InstrumentsSchema: InstrumentsSchema().index,
    MeasurementsSchema: MeasurementsSchema().index,
    PointsSchema: PointsSchema().index,
    StationsSchema: StationsSchema().index,
    TargetsSchema: TargetsSchema().index,
}

SCHEMA_COLUMNS = {
    "required": True,
    "unit": DistanceUnit.METER,
    "validator": IsNumeric(),
    "description": "A test column",
}

SCHEMA_EMPTY_COLUMNS = {
    column: value for column, value in SCHEMA_COLUMNS.items() if column != "required"
}


# ======================================================================================
# Tests of strict schema behavior
# ======================================================================================
@pytest.mark.parametrize("schema_class", SCHEMA_CLASSES_AND_INDEXES)
def test_validation_mode_after_initialization(schema_class):
    schema = schema_class()

    schema.validation_mode = "eager"

    # ----------------------------------------------------------------------------------

    assert schema.validation_mode is ValidationMode.EAGER


# ======================================================================================
@pytest.mark.parametrize("schema_class", SCHEMA_CLASSES_AND_INDEXES)
def test_adding_attributes(schema_class):
    schema = schema_class()

    # ----------------------------------------------------------------------------------

    with pytest.raises(AttributeError, match="does not allow modifying attribute"):
        schema.some_attribute = object()


# ======================================================================================
# Test setitem behavior for strict schema classes
# ======================================================================================
class TestSetItem:
    """Test setitem behavior for strict schema classes."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column", "value"),
        [
            (schema_class, column, value)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column, value in SCHEMA_COLUMNS.items()
        ],
    )
    def test_valid(self, schema_class, column, value):
        """Test that valid assignments are accepted and committed correctly."""

        schema = schema_class()
        dtypes_before = schema.dtypes.copy(deep=True)
        schema[column] = value

        # ------------------------------------------------------------------------------

        assert (schema[column] == value).all()
        assert schema.dtypes.equals(dtypes_before)

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column"),
        [
            (schema_class, column)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column in SCHEMA_EMPTY_COLUMNS
        ],
    )
    def test_empty(self, schema_class, column):
        """Test that setting a column to None is allowed."""
        schema = schema_class()
        schema[column] = None

        # ------------------------------------------------------------------------------

        assert schema[column].isnull().all()

    # ==================================================================================
    @pytest.mark.parametrize(
        "schema_class", [schema_class for schema_class in SCHEMA_CLASSES_AND_INDEXES]
    )
    def test_empty_required(self, schema_class):
        """Test that setting a None to 'required' column is not allowed."""
        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match="required"):
            schema["required"] = None

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column"),
        [
            (schema_class, column)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column in SCHEMA_COLUMNS
        ],
    )
    def test_invalid(self, schema_class, column):
        """Test that invalid assignments raise a ValueError."""

        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match=column):
            schema[column] = object()


# ======================================================================================
# Test loc behavior for strict schema classes
# ======================================================================================
class TestLoc:
    """Test label-based indexer behavior for strict schema classes."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column", "value"),
        [
            (schema_class, column, value)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column, value in SCHEMA_COLUMNS.items()
        ],
    )
    def test_valid(self, schema_class, column, value):
        """Test that valid column assignments are accepted and committed correctly."""

        schema = schema_class()
        dtypes_before = schema.dtypes.copy(deep=True)
        schema.loc[:, column] = value

        # ------------------------------------------------------------------------------

        assert (schema[column] == value).all()
        assert schema.dtypes.equals(dtypes_before)

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column"),
        [
            (schema_class, column)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column in SCHEMA_EMPTY_COLUMNS
        ],
    )
    def test_empty(self, schema_class, column):
        """Test that setting a column to None is allowed."""
        schema = schema_class()
        schema.loc[:, column] = None

        # ------------------------------------------------------------------------------

        assert schema[column].isnull().all()

    # ==================================================================================
    @pytest.mark.parametrize(
        "schema_class", [schema_class for schema_class in SCHEMA_CLASSES_AND_INDEXES]
    )
    def test_empty_required(self, schema_class):
        """Test that setting a None to 'required' column is not allowed."""
        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match="required"):
            schema.loc[:, "required"] = None

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column"),
        [
            (schema_class, column)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column in SCHEMA_COLUMNS
        ],
    )
    def test_invalid(self, schema_class, column):
        """Test that invalid assignments raise a ValueError."""

        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match=column):
            schema.loc[:, column] = object()

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column", "value"),
        [
            (schema_class, row, column, value)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column, value in SCHEMA_COLUMNS.items()
        ],
    )
    def test_single_item_valid(self, schema_class, row, column, value):
        """Test that loc supports valid single-item assignments."""

        schema = schema_class()
        schema.loc[row, column] = value

        # ------------------------------------------------------------------------------

        assert schema.loc[row, column] == value

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column"),
        [
            (schema_class, row, column)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column in SCHEMA_EMPTY_COLUMNS
        ],
    )
    def test_single_item_empty(self, schema_class, row, column):
        """Test that loc supports valid single-item assignments."""

        schema = schema_class()
        schema.loc[row, column] = None

        # ------------------------------------------------------------------------------

        assert pd.isna(schema.loc[row, column])

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row"),
        [
            (schema_class, row)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
        ],
    )
    def test_single_item_empty_required(self, schema_class, row):
        """Test that setting a None to 'required' column is not allowed."""

        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match="required"):
            schema.loc[row, "required"] = None

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column"),
        [
            (schema_class, row, column)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column in SCHEMA_COLUMNS
        ],
    )
    def test_single_item_invalid(self, schema_class, row, column):
        """Test that loc supports valid single-item assignments."""

        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match=column):
            schema.loc[row, column] = object()


# ======================================================================================
# Test iloc behavior for strict schema classes
# ======================================================================================
class TestIloc:
    """Test position-based indexer behavior for strict schema classes."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column", "value"),
        [
            (schema_class, column, value)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column, value in SCHEMA_COLUMNS.items()
        ],
    )
    def test_valid(self, schema_class, column, value):
        """Test that valid column assignments are accepted and committed correctly."""

        schema = schema_class()
        dtypes_before = schema.dtypes.copy(deep=True)
        column_position = schema.columns.get_loc(column)
        schema.iloc[:, column_position] = value

        # ------------------------------------------------------------------------------

        assert (schema[column] == value).all()
        assert schema.dtypes.equals(dtypes_before)

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column"),
        [
            (schema_class, column)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column in SCHEMA_EMPTY_COLUMNS
        ],
    )
    def test_empty(self, schema_class, column):
        """Test that setting a column to None is allowed."""
        schema = schema_class()
        column_position = schema.columns.get_loc(column)
        schema.iloc[:, column_position] = None

        # ------------------------------------------------------------------------------

        assert schema[column].isnull().all()

    # ==================================================================================
    @pytest.mark.parametrize(
        "schema_class", [schema_class for schema_class in SCHEMA_CLASSES_AND_INDEXES]
    )
    def test_empty_required(self, schema_class):
        """Test that setting a None to 'required' column is not allowed."""
        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match="required"):
            schema.loc[:, "required"] = None

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "column"),
        [
            (schema_class, column)
            for schema_class in SCHEMA_CLASSES_AND_INDEXES
            for column in SCHEMA_COLUMNS
        ],
    )
    def test_invalid(self, schema_class, column):
        """Test that invalid assignments raise a ValueError."""

        schema = schema_class()
        column_position = schema.columns.get_loc(column)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match=column):
            schema.iloc[:, column_position] = object()

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column", "value"),
        [
            (schema_class, row, column, value)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column, value in SCHEMA_COLUMNS.items()
        ],
    )
    def test_single_item_valid(self, schema_class, row, column, value):
        """Test that iloc supports valid single-item assignments."""

        schema = schema_class()
        dtypes_before = schema.dtypes.copy(deep=True)
        row_position = schema.index.get_loc(row)
        column_position = schema.columns.get_loc(column)
        schema.iloc[row_position, column_position] = value

        # ------------------------------------------------------------------------------

        assert schema.iloc[row_position, column_position] == value
        assert schema.dtypes.equals(dtypes_before)

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column"),
        [
            (schema_class, row, column)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column in SCHEMA_EMPTY_COLUMNS
        ],
    )
    def test_single_item_empty(self, schema_class, row, column):
        """Test that iloc supports valid single-item assignments."""

        schema = schema_class()
        row_position = schema.index.get_loc(row)
        column_position = schema.columns.get_loc(column)
        schema.iloc[row_position, column_position] = None

        # ------------------------------------------------------------------------------

        assert pd.isna(schema.iloc[row_position, column_position])

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row"),
        [
            (schema_class, row)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
        ],
    )
    def test_single_item_empty_required(self, schema_class, row):
        """Test that setting a None to 'required' column is not allowed."""

        schema = schema_class()
        row_position = schema.index.get_loc(row)
        column_position = schema.columns.get_loc("required")

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match="required"):
            schema.iloc[row_position, column_position] = None

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column"),
        [
            (schema_class, row, column)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column in SCHEMA_COLUMNS
        ],
    )
    def test_single_item_invalid(self, schema_class, row, column):
        """Test that iloc supports valid single-item assignments."""

        schema = schema_class()
        row_position = schema.index.get_loc(row)
        column_position = schema.columns.get_loc(column)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match=column):
            schema.iloc[row_position, column_position] = object()


# ======================================================================================
# Test at behavior for strict schema classes
# ======================================================================================
class TestAt:
    """Test scalar label-based indexer behavior for strict schema classes."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column", "value"),
        [
            (schema_class, row, column, value)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column, value in SCHEMA_COLUMNS.items()
        ],
    )
    def test_valid(self, schema_class, row, column, value):
        """Test that valid scalar assignments are accepted and committed correctly."""

        schema = schema_class()
        dtypes_before = schema.dtypes.copy(deep=True)
        schema.at[row, column] = value

        # ------------------------------------------------------------------------------

        assert schema.at[row, column] == value
        assert schema.dtypes.equals(dtypes_before)

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column"),
        [
            (schema_class, row, column)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column in SCHEMA_EMPTY_COLUMNS
        ],
    )
    def test_empty(self, schema_class, row, column):
        """Test that setting a scalar item to None is allowed."""
        schema = schema_class()
        schema.at[row, column] = None

        # ------------------------------------------------------------------------------

        assert pd.isna(schema.at[row, column])

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row"),
        [
            (schema_class, row)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
        ],
    )
    def test_empty_required(self, schema_class, row):
        """Test that setting a None to 'required' column is not allowed."""

        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match="required"):
            schema.at[row, "required"] = None

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column"),
        [
            (schema_class, row, column)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column in SCHEMA_COLUMNS
        ],
    )
    def test_invalid(self, schema_class, row, column):
        """Test that invalid scalar assignments raise a ValueError."""

        schema = schema_class()

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match=column):
            schema.at[row, column] = object()


# ======================================================================================
# Test iat behavior for strict schema classes
# ======================================================================================
class TestIat:
    """Test scalar position-based indexer behavior for strict schema classes."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column", "value"),
        [
            (schema_class, row, column, value)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column, value in SCHEMA_COLUMNS.items()
        ],
    )
    def test_valid(self, schema_class, row, column, value):
        """Test that valid scalar assignments are accepted and committed correctly."""

        schema = schema_class()
        dtypes_before = schema.dtypes.copy(deep=True)
        column_position = schema.columns.get_loc(column)
        row_position = schema.index.get_loc(row)
        schema.iat[row_position, column_position] = value

        # ------------------------------------------------------------------------------

        assert schema.iat[row_position, column_position] == value
        assert schema.dtypes.equals(dtypes_before)

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column"),
        [
            (schema_class, row, column)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column in SCHEMA_EMPTY_COLUMNS
        ],
    )
    def test_empty(self, schema_class, row, column):
        """Test that setting a scalar item to None is allowed."""
        schema = schema_class()
        column_position = schema.columns.get_loc(column)
        row_position = schema.index.get_loc(row)
        schema.iat[row_position, column_position] = None

        # ------------------------------------------------------------------------------

        assert pd.isna(schema.iat[row_position, column_position])

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row"),
        [
            (schema_class, row)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
        ],
    )
    def test_empty_required(self, schema_class, row):
        """Test that setting a None to 'required' column is not allowed."""

        schema = schema_class()
        row_position = schema.index.get_loc(row)
        column_position = schema.columns.get_loc("required")

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match="required"):
            schema.iat[row_position, column_position] = None

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "row", "column"),
        [
            (schema_class, row, column)
            for schema_class, index in SCHEMA_CLASSES_AND_INDEXES.items()
            for row in index
            for column in SCHEMA_COLUMNS
        ],
    )
    def test_invalid(self, schema_class, row, column):
        """Test that invalid scalar assignments raise a ValueError."""

        schema = schema_class()
        column_position = schema.columns.get_loc(column)
        row_position = schema.index.get_loc(row)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValueError, match=column):
            schema.iat[row_position, column_position] = object()
