# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from abc import abstractmethod
from typing import Any, Generic, Self, cast

from pysurveyor.abc.pysurveyor_abc import ABSTRACT_CLASS_ATTRIBUTE, PySurveyorABC
from pysurveyor.io.accessors.loader_accessor import LoaderAccessor
from pysurveyor.typing.vars import LoaderType


# ======================================================================================
# Abstract base class for PySurveyor I/O objects
# ======================================================================================
class PySurveyorIO(PySurveyorABC, Generic[LoaderType]):
    """Base class for PySurveyor objects that support source loading.

    ``PySurveyorIO`` requires concrete subclasses to provide a
    :class:`LoaderAccessor` through the ``load`` class attribute. The accessor
    exposes a loader bound to the owning domain class, allowing loading logic.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    load = cast(LoaderAccessor[LoaderType], ABSTRACT_CLASS_ATTRIBUTE)

    # ==================================================================================
    @abstractmethod
    def copy(self, deep: bool = True) -> Self:
        """Return a copy of the table.

        Parameters
        ----------
        deep : bool, default=True
            Whether to copy the underlying data deeply.

        Returns
        -------
        Self
            A copy of the table.
        """
        pass

    # ==================================================================================
    @abstractmethod
    def _copy_kwargs(self) -> dict[str, Any]:
        """Return extra arguments for a domain-specific copy constructor.

        Returns
        -------
        dict[str, Any]
            Keyword arguments passed to the concrete table constructor when
            creating a copy.
        """
        pass
