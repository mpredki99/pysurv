# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests for I/O configuration objects."""

from pysurveyor.io.config import CSVConfig


# ======================================================================================
# Test CSV file configuration class
# ======================================================================================
class TestCSVConfig:
    """Test CSV configuration construction and access."""

    # ==================================================================================
    def test_csv_config_kwargs(self):
        """Expose configured CSV options through ``kwargs``."""
        skiprows = lambda index: index % 2 == 0

        config = CSVConfig(
            delimiter=";",
            decimal=",",
            skiprows=skiprows,
            nrows=10,
            comment="#",
        )

        # ------------------------------------------------------------------------------

        assert config.kwargs == {
            "delimiter": ";",
            "decimal": ",",
            "skiprows": skiprows,
            "nrows": 10,
            "comment": "#",
        }

    # ==================================================================================
    def test_kwargs_is_fresh(self):
        """Return a fresh options dictionary for each access."""
        config = CSVConfig()

        first_kwargs = config.kwargs
        second_kwargs = config.kwargs

        # ------------------------------------------------------------------------------

        assert first_kwargs == second_kwargs
        assert first_kwargs is not second_kwargs

    # ==================================================================================
    def test_update_config_updates_kwargs(self):
        """Keep ``kwargs`` synchronized with updated configuration fields."""
        config = CSVConfig(delimiter=";")

        assert config.delimiter == ";"
        assert config.kwargs["delimiter"] == ";"

        # ------------------------------------------------------------------------------

        config.delimiter = ","

        assert config.delimiter == ","
        assert config.kwargs["delimiter"] == ","
