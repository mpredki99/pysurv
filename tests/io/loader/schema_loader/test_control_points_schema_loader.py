"""Tests for the ControlPointsSchema loader."""

import pytest

from pysurveyor.schema import ControlPointsSchema
from pysurveyor.schema.model import models


# ======================================================================================
# Test control points schema loader from CSV files
# ======================================================================================
class TestFromCSV:

    # ==================================================================================
    # Control points schema data frame fixture
    # ==================================================================================
    control_points_schema_df = pytest.fixture(
        lambda self, test_data_frame: test_data_frame("control_points_schema.csv")
    )

    # ==================================================================================
    def test_read_invalid_path(self):
        with pytest.raises(FileNotFoundError):
            ControlPointsSchema.load.from_csv("non_existent_file.csv")

    # ==================================================================================
    def test_read_valid_csv(self, control_points_schema_df, csv_file):
        # Write the DataFrame to a temporary CSV file ----------------------------------
        df = control_points_schema_df.valid
        csv_filepath = csv_file(df)

        # Load the schema from the CSV file --------------------------------------------
        schema = ControlPointsSchema.load.from_csv(csv_filepath)

        # ------------------------------------------------------------------------------
        assert schema.index.equals(models.control_points.index)
        assert schema.columns.equals(models.control_points.columns)
        assert schema.dtypes.equals(models.control_points.dtypes)

    # ==================================================================================
    def test_missing_index_column(self, control_points_schema_df, csv_file):
        # Write the DataFrame to a temporary CSV file ----------------------------------
        df = control_points_schema_df.missing_column("field")
        csv_filepath = csv_file(df)

        # Load the schema from the CSV file --------------------------------------------
        with pytest.raises(ValueError, match=f"Invalid index"):
            ControlPointsSchema.load.from_csv(csv_filepath)

    # ==================================================================================
    @pytest.mark.parametrize(
        "missing_column",
        [
            "required",
            "unit",
            "validator",
            "description",
        ],
    )
    def test_missing_column(self, control_points_schema_df, missing_column, csv_file):
        # Write the DataFrame to a temporary CSV file ----------------------------------
        df = control_points_schema_df.missing_column(missing_column)
        csv_filepath = csv_file(df)

        # Load the schema from the CSV file --------------------------------------------
        with pytest.raises(ValueError, match=f"{missing_column}"):
            ControlPointsSchema.load.from_csv(csv_filepath)
