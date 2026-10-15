# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Generic

from pysurveyor.abc.pysurveyor_abc import PySurveyorABC
from pysurveyor.typing.vars import ConfigType


# ======================================================================================
# Abstract base class for data sources
# ======================================================================================
class Source(PySurveyorABC):
    """Base class for objects describing where data comes from.

    A source stores location and format-specific configuration.
    """

    __pysurveyor_abstract__ = True


# ======================================================================================
# File source
# ======================================================================================
@dataclass
class FileSource(Source, Generic[ConfigType]):
    """Source describing a file path and its configuration."""

    __pysurveyor_abstract__ = True

    path: Path
    config: ConfigType
