"""Tests for the common table-loader orchestration."""

from pathlib import Path
from typing import Any

from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.loader.abc.table_loader import TableLoader
from pysurveyor.io.source.file_sources import CSVSource


# ======================================================================================
# Implementation of necessary classes for testing purposes
# ======================================================================================
class DomainObject:
    """Minimal domain object for testing purposes."""


# ======================================================================================
class ConcreteTableLoader(TableLoader[DomainObject]):  # type: ignore
    """Test implementation that records the arguments passed to ``from_source``."""

    # ==================================================================================
    def __init__(self) -> None:
        super().__init__(DomainObject)
        self.received_source = None
        self.received_kwargs = None

    # ==================================================================================
    def from_source(self, source: CSVSource, **kwargs: Any) -> DomainObject:
        self.received_source = source
        self.received_kwargs = kwargs
        return DomainObject()

    # ==================================================================================
    # Minimal implementation of abstract method
    # ==================================================================================
    def from_csv(self, *args: Any, **kwargs: Any) -> DomainObject:
        raise NotImplementedError


loader = ConcreteTableLoader()


# ======================================================================================
# Test functions
# ======================================================================================
def test_passing_existing_source(tmp_path: Path) -> None:
    # ----------------------------------------------------------------------------------
    # Create a CSV file and a corresponding CSVSource
    # ----------------------------------------------------------------------------------
    path = tmp_path / "data.csv"
    path.write_text("name,value\nexample,1\n", encoding="utf-8")

    source = CSVSource(path=path, config=CSVConfig())

    # ----------------------------------------------------------------------------------
    loader(source, option="value")

    assert loader.received_source is source
    assert loader.received_kwargs == {"option": "value"}


# ======================================================================================
def test_call_resolves_source(tmp_path: Path) -> None:
    # ----------------------------------------------------------------------------------
    # Create a CSV config and a CSV file for testing
    # ----------------------------------------------------------------------------------
    path = tmp_path / "data.csv"
    path.write_text("name,value\nexample,1\n", encoding="utf-8")

    config = CSVConfig(delimiter=";")

    # ----------------------------------------------------------------------------------
    loader(path, config=config, option="value")

    assert isinstance(loader.received_source, CSVSource)
    assert loader.received_source.path == path.resolve()
    assert loader.received_source.config is config
    assert loader.received_kwargs == {"option": "value"}
