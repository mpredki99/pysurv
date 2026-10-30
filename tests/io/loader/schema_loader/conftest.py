from pathlib import Path

import pandas as pd
import pytest

# ======================================================================================
# Directory containing CSV files for schema tests
# ======================================================================================
CSV_DIRECTORY = Path(__file__).parents[3] / "fixtures" / "csv" / "schema"


# ======================================================================================
# Test DataFrame fixture
# ======================================================================================
@pytest.fixture
def test_data_frame():

    # ----------------------------------------------------------------------------------
    # Create a class to handle test DataFrames for tests
    # ----------------------------------------------------------------------------------
    class TestDataFrame:

        # ------------------------------------------------------------------------------
        def __init__(self, filename: str):
            self.filename = filename

        # ------------------------------------------------------------------------------
        @property
        def valid(self):
            """Return a valid DataFrame for the given filename."""
            filepath = CSV_DIRECTORY / self.filename
            return pd.read_csv(filepath)

        # ------------------------------------------------------------------------------
        def missing_column(self, column: str) -> pd.DataFrame:
            """Return a DataFrame with a missing column."""
            return self.valid.drop(columns=[column])

        # ------------------------------------------------------------------------------
        def missing_row(self, row_index: int) -> pd.DataFrame:
            """Return a DataFrame with a missing row."""
            return self.valid.drop(index=[row_index])

        # ------------------------------------------------------------------------------
        def add_column(self, column: str, value=None) -> pd.DataFrame:
            """Return a DataFrame with an additional column."""
            result = self.valid.copy()
            result[column] = value
            return result

        # ------------------------------------------------------------------------------
        def add_row(self, row: dict) -> pd.DataFrame:
            """Return a DataFrame with an additional row."""
            return pd.concat(
                [self.valid.copy(), pd.DataFrame([row])], ignore_index=True
            )

        # ------------------------------------------------------------------------------
        def replace_value(self, row_index: int, column: str, value) -> pd.DataFrame:
            """Return a DataFrame with a replaced value."""
            result = self.valid.copy()
            result.at[row_index, column] = value
            return result

        # ------------------------------------------------------------------------------
        def replace_with_empty(self, row_index: int, column: str) -> pd.DataFrame:
            """Return a DataFrame with a replaced empty value."""
            return self.replace_value(row_index, column, pd.NA)

    # ----------------------------------------------------------------------------------
    return TestDataFrame


# ======================================================================================
# Temporary CSV file fixture
# ======================================================================================
@pytest.fixture
def csv_file(tmp_path):
    """Write a DataFrame to a temporary CSV file."""

    def _create(df, filename="test.csv"):
        filepath = tmp_path / filename
        df.to_csv(filepath, index=False)
        return filepath

    return _create
