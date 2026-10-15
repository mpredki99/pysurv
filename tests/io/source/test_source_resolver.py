"""Tests for source resolution."""

from os import PathLike
from pathlib import Path

import pytest

from pysurveyor.io.config.abc.config import FileConfig
from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.source.file_sources import CSVSource
from pysurveyor.io.source.resolver import resolve_source


# ======================================================================================
# Function to create a sample CSV file for testing
# ======================================================================================
def create_csv_file(directory: Path, name: str = "data.csv") -> Path:
    """Create a sample CSV file."""
    path = directory / name
    path.write_text("name,value\nexample,1\n", encoding="utf-8")
    return path


# ======================================================================================
# Test functions for CSVSource resolving
# ======================================================================================
def test_resolve_string(tmp_path: Path) -> None:
    path = create_csv_file(tmp_path)

    source = resolve_source(str(path))

    assert isinstance(source, CSVSource)
    assert isinstance(source.config, CSVConfig)
    assert source.path == path.resolve()


# ======================================================================================
def test_resolve_path(tmp_path: Path) -> None:
    path = create_csv_file(tmp_path)

    source = resolve_source(path)

    assert isinstance(source, CSVSource)
    assert isinstance(source.config, CSVConfig)
    assert source.path == path.resolve()


# ======================================================================================
def test_resolve_pathlike(tmp_path: Path) -> None:
    # ----------------------------------------------------------------------------------
    # Create a path-like object
    # ----------------------------------------------------------------------------------
    class StringPathLike(PathLike[str]):
        def __init__(self, path: str) -> None:
            self.path = path

        def __fspath__(self) -> str:
            return self.path

    # ----------------------------------------------------------------------------------

    path = create_csv_file(tmp_path)

    source_input = StringPathLike(str(path))

    source = resolve_source(source_input)

    assert isinstance(source, CSVSource)
    assert isinstance(source.config, CSVConfig)
    assert source.path == path.resolve()


# ======================================================================================
def test_resolve_with_config(tmp_path: Path) -> None:
    path = create_csv_file(tmp_path)
    config = CSVConfig(delimiter=";", decimal=",", nrows=1)

    source = resolve_source(path, config=config)

    assert isinstance(source, CSVSource)
    assert source.config is config


# ======================================================================================
def test_resolve_uppercase_suffix(tmp_path: Path) -> None:
    path = create_csv_file(tmp_path, "data.CSV")

    source = resolve_source(path)

    assert isinstance(source, CSVSource)
    assert source.path == path.resolve()


# ======================================================================================
def test_resolve_relative_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = create_csv_file(tmp_path)
    monkeypatch.chdir(tmp_path)

    source = resolve_source(Path(path.name))

    assert source.path == path.resolve()  # type: ignore


# ======================================================================================
def test_reject_unsupported_input() -> None:
    with pytest.raises(TypeError, match=r"Could not resolve source from input"):
        resolve_source(123)


# ======================================================================================
def test_reject_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match=r"does not exist"):
        resolve_source(tmp_path / "missing.csv")


# ======================================================================================
def test_reject_unsupported_suffix(tmp_path: Path) -> None:
    path = tmp_path / "data.txt"
    path.write_text("data", encoding="utf-8")

    with pytest.raises(TypeError, match=r"Unsupported file suffix"):
        resolve_source(path)


# ======================================================================================
def test_reject_incompatible_config(tmp_path: Path) -> None:
    # ----------------------------------------------------------------------------------
    # Create custom file configuration class that is incompatible with CSVSource
    # ----------------------------------------------------------------------------------
    class OtherFileConfig(FileConfig):
        """Configuration incompatible with CSV sources."""

        option: str = "value"

    path = create_csv_file(tmp_path)

    with pytest.raises(TypeError, match=r"CSV sources require CSVConfig"):
        resolve_source(path, config=OtherFileConfig())
