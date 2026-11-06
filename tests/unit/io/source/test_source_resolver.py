# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests for source resolution."""

from os import PathLike
from pathlib import Path

import pytest

from pysurveyor.io.config.abc.config import FileConfig
from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.source.file_sources import CSVSource
from pysurveyor.io.source.resolver import resolve_source


# ======================================================================================
# Implementation of a PathLike class for testing purposes
# ======================================================================================
class StringPathLike(PathLike[str]):
    def __init__(self, path):
        self.path = str(path)

    def __fspath__(self):
        return self.path


# ======================================================================================
# General tests for resolve_source function
# ======================================================================================
def test_unsupported_input():
    """Test that resolve_source raises TypeError for unsupported input types."""
    input = 123

    # ----------------------------------------------------------------------------------

    with pytest.raises(TypeError, match=str(input)):
        resolve_source(input)


# ======================================================================================
def test_missing_file(tmp_path):
    input = tmp_path / "missing.csv"

    # ----------------------------------------------------------------------------------

    with pytest.raises(FileNotFoundError, match=str(input.name)):
        resolve_source(input)


# ======================================================================================
def test_unsupported_suffix(tmp_text_file):

    input = tmp_text_file("data.unsupported", data="unsupported suffix")

    # ----------------------------------------------------------------------------------

    with pytest.raises(TypeError, match=str(input.suffix)):
        resolve_source(input)


# ======================================================================================
# Test functions for CSVSource resolving
# ======================================================================================
class TestResolveCSVSource:
    """Tests for resolving CSV sources."""

    @pytest.mark.parametrize(
        "path_factory",
        [
            pytest.param(str, id="string"),
            pytest.param(Path, id="path"),
            pytest.param(StringPathLike, id="pathlike"),
        ],
    )
    def test_resolve(self, path_factory, tmp_text_file):
        """Test that resolve_source correctly resolves CSV files."""
        filepath = tmp_text_file("test.csv")
        source = resolve_source(path_factory(filepath))

        # ------------------------------------------------------------------------------

        assert isinstance(source, CSVSource)
        assert isinstance(source.config, CSVConfig)
        assert source.path == filepath.resolve()

    # ==================================================================================
    def test_resolve_with_config(self, tmp_text_file):
        """Test that resolve_source uses a provided CSVConfig."""

        filepath = tmp_text_file("test.csv")
        config = CSVConfig(delimiter=";", decimal=",", nrows=1)

        source = resolve_source(filepath, config=config)

        # ------------------------------------------------------------------------------

        assert isinstance(source, CSVSource)
        assert source.config is config

    # ==================================================================================
    def test_uppercase_suffix(self, tmp_text_file):
        """Test that resolve_source can handle uppercase file suffixes."""

        filepath = tmp_text_file("test.CSV")
        source = resolve_source(filepath)

        # ------------------------------------------------------------------------------

        assert isinstance(source, CSVSource)
        assert source.path == filepath.resolve()

    # ==================================================================================
    def test_relative_path(self, tmp_text_file, monkeypatch: pytest.MonkeyPatch):
        """Test that resolve_source can handle relative paths."""

        filepath = tmp_text_file("test.csv")
        monkeypatch.chdir(filepath.parent)

        source = resolve_source(Path(filepath.name))

        # ------------------------------------------------------------------------------

        assert source.path == filepath.resolve()  # type: ignore

    # ==================================================================================
    def test_incompatible_config(self, tmp_text_file):
        """Test that resolve_source raises TypeError for incompatible config types."""

        # ------------------------------------------------------------------------------
        # Custom file configuration class that is incompatible with CSVSource
        # ------------------------------------------------------------------------------
        class OtherFileConfig(FileConfig):
            option: str = "value"

        # ------------------------------------------------------------------------------

        filepath = tmp_text_file("test.csv")

        # ------------------------------------------------------------------------------

        with pytest.raises(TypeError, match=r"CSV sources require CSVConfig"):
            resolve_source(filepath, config=OtherFileConfig())
