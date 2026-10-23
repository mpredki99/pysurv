# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from pysurveyor.io.config.abc.config import FileConfig


# ======================================================================================
# CSV file configuration
# ======================================================================================
@dataclass
class CSVConfig(FileConfig):
    """Configuration options for reading CSV files."""

    delimiter: str | None = None
    decimal: str = "."
    skiprows: int | Callable | None = None
    nrows: int | None = None
    comment: str | None = "#"
