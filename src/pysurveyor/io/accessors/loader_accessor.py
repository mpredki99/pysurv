# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Generic, cast

from pysurveyor.typing.vars import LoaderType

if TYPE_CHECKING:
    from pysurveyor.io.abc.pysurveyor_io import PySurveyorIO


# ======================================================================================
# Class descriptor for constructor-bound loaders
# ======================================================================================
class LoaderAccessor(Generic[LoaderType]):
    """Expose a loader bound to its owning domain class.

    Each owning class caches its own loader in its class dictionary.
    Subclasses therefore receive loaders bound to their own class instead
    of reusing a loader bound to a parent class.

    Attributes
    ----------
    _constructor : type[LoaderType]
        Loader class used to create the bound loader instance.
    _attr : str
        Name of the class attribute used to cache the loader instance.
    """

    # ==================================================================================
    def __init__(self, constructor: type[LoaderType]) -> None:
        """Initialize the accessor with the loader constructor."""
        self._constructor = constructor
        self._attr = ""

    # ==================================================================================
    def __set_name__(self, owner: type[Any], name: str) -> None:
        """Register the class attribute used to cache the loader."""
        self._attr = f"_loader_of_{name}"

    # ==================================================================================
    def __get__(
        self,
        obj: PySurveyorIO[Any] | None,
        owner: type[PySurveyorIO[Any]],
    ) -> LoaderType:
        """Return the loader bound to the owning class.

        The loader is created on first access and reused for subsequent
        class or instance accesses.

        Parameters
        ----------
        obj : PySurveyorIO | None
            Instance through which the descriptor is accessed.
        owner : type[PySurveyorIO]
            Class that owns the loader binding.

        Returns
        -------
        LoaderType
            Loader instance constructed with ``owner``.
        """
        if self._attr not in owner.__dict__:
            setattr(owner, self._attr, self._constructor(owner))

        return cast(LoaderType, owner.__dict__[self._attr])
