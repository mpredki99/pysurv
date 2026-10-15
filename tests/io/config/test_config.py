"""Tests for I/O configuration objects."""

from pysurveyor.io.config.file_configs import CSVConfig


# ======================================================================================
# Test `kwargs` property
# ======================================================================================
def test_csv_config_kwargs() -> None:
    skiprows = lambda index: index % 2 == 0

    config = CSVConfig(
        delimiter=";",
        decimal=",",
        skiprows=skiprows,
        nrows=10,
    )

    assert config.kwargs == {
        "delimiter": ";",
        "decimal": ",",
        "skiprows": skiprows,
        "nrows": 10,
    }


# ======================================================================================
def test_kwargs_is_fresh() -> None:
    config = CSVConfig()

    first_kwargs = config.kwargs
    second_kwargs = config.kwargs

    assert first_kwargs == second_kwargs
    assert first_kwargs is not second_kwargs


# ======================================================================================
def test_update_config_updates_kwargs() -> None:
    config = CSVConfig(delimiter=";")

    assert config.delimiter == ";"
    assert config.kwargs["delimiter"] == ";"

    config.delimiter = ","

    assert config.delimiter == ","
    assert config.kwargs["delimiter"] == ","
