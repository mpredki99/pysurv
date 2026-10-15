# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from abc import abstractmethod
from dataclasses import dataclass, fields
from typing import Any

from pysurveyor.abc.pysurveyor_abc import PySurveyorABC


# ======================================================================================
# Abstract base class for I/O configurations
# ======================================================================================
class IOConfig(PySurveyorABC):
    """Base class for configurations used by PySurveyor I/O operations.

    ``IOConfig`` defines the interface for converting configuration state into
    keyword arguments accepted by a concrete backend reader.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    @property
    @abstractmethod
    def kwargs(self) -> dict[str, Any]:
        """Return backend keyword arguments for the I/O operation."""
        pass


# ======================================================================================
# File configuration classes
# ======================================================================================
@dataclass
class FileConfig(IOConfig):
    """Base configuration for file-based I/O operations."""

    __pysurveyor_abstract__ = True

    # ==================================================================================
    @property
    def kwargs(self) -> dict[str, Any]:
        """Return a fresh mapping of configuration fields for the backend.

        Values are returned unchanged, including callable values and other
        backend-specific objects. A shallow snapshot prevents later
        configuration changes from modifying an already-created kwargs mapping.
        """
        return {field.name: getattr(self, field.name) for field in fields(self)}
