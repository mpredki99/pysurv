"""Tests of public API for loading schema definitions."""

from pathlib import Path

import pytest

from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.source import CSVSource
from pysurveyor.schema import (
    ControlPointsSchema,
    InstrumentsSchema,
    MeasurementsSchema,
    PointsSchema,
    StationsSchema,
    TargetsSchema,
    models,
)

# ======================================================================================
# Path to the directory containing the schema CSV files
# ======================================================================================
SCHEMA_CSV_DIR = Path(__file__).resolve().parents[1] / "data" / "csv" / "schema"

# ======================================================================================
# Test data for schema loading tests
# ======================================================================================
SCHEMA_TEST_DATA = {
    ControlPointsSchema: {
        "model": models.control_points,
        "csv_file": SCHEMA_CSV_DIR / "control_points_schema.csv",
    },
    InstrumentsSchema: {
        "model": models.instruments,
        "csv_file": SCHEMA_CSV_DIR / "instruments_schema.csv",
    },
    MeasurementsSchema: {
        "model": models.measurements,
        "csv_file": SCHEMA_CSV_DIR / "measurements_schema.csv",
    },
    PointsSchema: {
        "model": models.points,
        "csv_file": SCHEMA_CSV_DIR / "points_schema.csv",
    },
    StationsSchema: {
        "model": models.stations,
        "csv_file": SCHEMA_CSV_DIR / "stations_schema.csv",
    },
    TargetsSchema: {
        "model": models.targets,
        "csv_file": SCHEMA_CSV_DIR / "targets_schema.csv",
    },
}


# ======================================================================================
# Test classes for schema loading tests
# ======================================================================================


class TestLoad:
    """Verify the top-level ``Schema.load`` public entry point.

    These tests exercise automatic source resolution from user-provided input.
    Successful loads must produce the expected concrete schema structure, while
    unsupported input types must fail before loading begins.
    """

    # ==================================================================================
    @pytest.mark.parametrize(
        "schema_class",
        [schema_class for schema_class in SCHEMA_TEST_DATA.keys()],
    )
    def test_load_schema_from_invalid_type(self, schema_class):
        """Reject input types that cannot be resolved as schema sources."""
        # ------------------------------------------------------------------------------
        with pytest.raises(TypeError):
            schema_class.load(None)

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "model", "csv_file"),
        [
            (schema_class, data["model"], data["csv_file"])
            for schema_class, data in SCHEMA_TEST_DATA.items()
        ],
    )
    def test_load_schema_from_csv(self, schema_class, model, csv_file):
        """Load each concrete schema from its real CSV fixture."""
        # ------------------------------------------------------------------------------
        # Assert that test fixture CSV file exists
        # ------------------------------------------------------------------------------
        assert csv_file.is_file(), f"Missing test fixture: {csv_file.filename}"

        # ------------------------------------------------------------------------------
        loaded_schema = schema_class.load(csv_file)

        # ------------------------------------------------------------------------------
        assert loaded_schema.index.equals(model.index)
        assert loaded_schema.columns.equals(model.columns)
        assert loaded_schema.dtypes.equals(model.dtypes)


# ======================================================================================
class TestLoadFromSource:
    """Verify loading schemas from explicit ``Source`` objects.

    These tests cover the public boundary where callers construct a source
    object themselves and pass it to the schema loader.
    """

    # ==================================================================================
    @pytest.mark.parametrize(
        "schema_class",
        [schema_class for schema_class in SCHEMA_TEST_DATA.keys()],
    )
    def test_load_schema_from_invalid_type(self, schema_class):
        """Reject values that are not valid source objects."""
        # ------------------------------------------------------------------------------
        with pytest.raises(TypeError):
            schema_class.load.from_source(None)

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "model", "csv_file"),
        [
            (schema_class, data["model"], data["csv_file"])
            for schema_class, data in SCHEMA_TEST_DATA.items()
        ],
    )
    def test_load_schema_from_source(self, schema_class, model, csv_file):
        """Load each concrete schema through an explicit CSV source."""
        # ------------------------------------------------------------------------------
        # Assert that test fixture CSV file exists
        # ------------------------------------------------------------------------------
        assert csv_file.is_file(), f"Missing test fixture: {csv_file.name}"

        # ------------------------------------------------------------------------------
        source = CSVSource(path=csv_file)
        loaded_schema = schema_class.load.from_source(source)

        # ------------------------------------------------------------------------------
        assert loaded_schema.index.equals(model.index)
        assert loaded_schema.columns.equals(model.columns)
        assert loaded_schema.dtypes.equals(model.dtypes)


# ======================================================================================
class TestLoadFromCSV:
    """Verify the explicit ``Schema.load.from_csv`` public entry point.

    These tests cover direct path-based loading, including the public error
    raised for a missing file and the structure preserved for valid fixtures.
    """

    # ==================================================================================
    @pytest.mark.parametrize(
        "schema_class",
        [schema_class for schema_class in SCHEMA_TEST_DATA.keys()],
    )
    def test_load_schema_from_invalid_path(self, schema_class):
        """Raise ``FileNotFoundError`` when the CSV path does not exist."""
        # ------------------------------------------------------------------------------
        with pytest.raises(FileNotFoundError):
            schema_class.load.from_csv("invalid_path.csv")

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "model", "csv_file"),
        [
            (schema_class, data["model"], data["csv_file"])
            for schema_class, data in SCHEMA_TEST_DATA.items()
        ],
    )
    def test_load_schema_from_valid_path(self, schema_class, model, csv_file):
        """Load each concrete schema from a valid CSV path."""
        # ------------------------------------------------------------------------------
        # Assert that test fixture CSV file exists
        # ------------------------------------------------------------------------------
        assert csv_file.is_file(), f"Missing test fixture: {csv_file.name}"

        # ------------------------------------------------------------------------------
        loaded_schema = schema_class.load.from_csv(csv_file)

        # ------------------------------------------------------------------------------
        assert loaded_schema.index.equals(model.index)
        assert loaded_schema.columns.equals(model.columns)
        assert loaded_schema.dtypes.equals(model.dtypes)
        assert loaded_schema.validation_mode == "lazy"

    # ==================================================================================
    @pytest.mark.parametrize(
        ("schema_class", "model", "csv_file"),
        [
            (schema_class, data["model"], data["csv_file"])
            for schema_class, data in SCHEMA_TEST_DATA.items()
        ],
    )
    def test_load_schema_with_kwargs(self, schema_class, model, csv_file):
        """Load each concrete schema from a valid CSV path."""
        # ------------------------------------------------------------------------------
        # Assert that test fixture CSV file exists
        # ------------------------------------------------------------------------------
        assert csv_file.is_file(), f"Missing test fixture: {csv_file.name}"

        # ------------------------------------------------------------------------------
        config = CSVConfig(skiprows=lambda x: False)
        loaded_schema = schema_class.load.from_csv(
            csv_file, config=config, validation_mode="disabled"
        )

        # ------------------------------------------------------------------------------
        assert loaded_schema.index.equals(model.index)
        assert loaded_schema.columns.equals(model.columns)
        assert loaded_schema.dtypes.equals(model.dtypes)
        assert loaded_schema.validation_mode == "disabled"
