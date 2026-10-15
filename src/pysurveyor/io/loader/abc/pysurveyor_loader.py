# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from __future__ import annotations

from abc import abstractmethod
from typing import Any, Generic

from pysurveyor.abc.pysurveyor_abc import PySurveyorABC
from pysurveyor.typing.vars import DomainClassType


# ======================================================================================
# Abstract base class for loaders in PySurveyor
# ======================================================================================
class PySurveyorLoader(PySurveyorABC, Generic[DomainClassType]):
    """Base class for PySurveyor data loaders.

    ``PySurveyorLoader`` defines the interface for loading source data into
    domain objects.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    @abstractmethod
    def __call__(self, *args: Any, **kwargs: Any) -> DomainClassType:
        """Load data, resolving the source into source object when needed."""
        pass

    # ==================================================================================
    @abstractmethod
    def from_source(self, *args: Any, **kwargs: Any) -> DomainClassType:
        """Resolve a data source and build a domain object."""
        pass

    # ==================================================================================
    @abstractmethod
    def from_csv(self, *args: Any, **kwargs: Any) -> DomainClassType:
        """Import data from CSV files."""
        pass
