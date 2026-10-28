# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from contextlib import AbstractContextManager
from contextvars import ContextVar, Token
from os import PathLike
from pathlib import Path
from types import TracebackType
from typing import Callable, ClassVar, Hashable

from pysurveyor.backend.backend import DataFrame
from pysurveyor.io.cache.cache_keys import FileCacheKey


# ======================================================================================
# Loader cache to avoid repeated reads of the same source data
# ======================================================================================
class LoaderCache:
    """Cache loaded data within an execution context.

    Caching is disabled by default. Set LoaderCache.enabled = True to enable caching for
    the shared cache, or use session() for scoped caching that is automatically cleared
    on exit.
    """

    # Public cache switch --------------------------------------------------------------
    enabled: ClassVar[bool] = False
    # ----------------------------------------------------------------------------------

    _shared_cache: ClassVar[LoaderCache | None] = None
    _active_cache: ClassVar[ContextVar[LoaderCache | None]] = ContextVar(
        "current_loader_cache", default=None
    )

    # ==================================================================================
    def __init__(self, enabled: bool | None = None) -> None:
        """Create an empty cache."""
        self._enabled = enabled
        self._cached_data: dict[Hashable, DataFrame] = {}

    # ==================================================================================
    @property
    def is_enabled(self) -> bool:
        """Return whether caching is enabled for this instance."""
        if self._enabled is not None:
            return self._enabled
        return self.__class__.enabled

    # ==================================================================================
    # Class methods
    # ==================================================================================
    @classmethod
    def shared(cls) -> LoaderCache:
        """Return the shared instance of a cache."""
        if cls._shared_cache is None:
            cls._shared_cache = cls()
        return cls._shared_cache

    # ==================================================================================
    @classmethod
    def current(cls) -> LoaderCache:
        """Return the active cache, or fallback to the shared cache."""
        cache = cls._active_cache.get()
        return cache if cache is not None else cls.shared()

    # ==================================================================================
    @classmethod
    def session(cls) -> _CacheSession:
        """Create an isolated, automatically cleared cache session."""
        return _CacheSession(cls(enabled=True))

    # ==================================================================================
    # Public methods
    # ==================================================================================
    def get_data(
        self,
        key: Hashable,
        loader: Callable[[], DataFrame],
        *,
        copy: bool = True,
    ) -> DataFrame:
        """Load a data frame from the cache or use the loader to read it.

        Parameters:
        ----------
        key : Hashable
            Unique key to identify the data frame in the cache.
        loader : Callable[[], DataFrame]
            Function to load the data frame if it is not already cached.
        copy : bool, optional
            Return a copy of the cached frame when ``True`` (default).

        Returns:
        -------
        DataFrame
            The loaded data frame, either from the cache or freshly loaded.
        """
        if not self.is_enabled:
            return loader()

        # ------------------------------------------------------------------------------
        if key not in self._cached_data:
            self._cached_data[key] = loader()

        # ------------------------------------------------------------------------------
        data = self._cached_data[key]

        # ------------------------------------------------------------------------------
        return data.copy() if copy else data

    # ==================================================================================
    def clear(self) -> None:
        """Remove all cached data."""
        self._cached_data.clear()

    # ==================================================================================
    def remove(self, key: Hashable) -> None:
        """Remove one cache entry by key."""
        self._cached_data.pop(key, None)

    # ==================================================================================
    def remove_file(self, path: PathLike[str] | str) -> None:
        """Remove all file-cache entries associated with ``path``."""
        resolved_path = str(Path(path).expanduser().resolve())

        # ------------------------------------------------------------------------------
        self._cached_data = {
            key: data
            for key, data in self._cached_data.items()
            if not (
                isinstance(key, FileCacheKey) and key.resolved_path == resolved_path
            )
        }


# ======================================================================================
# Session context manager for a temporary loader cache
# ======================================================================================
class _CacheSession(AbstractContextManager[LoaderCache]):
    """Context manager that temporarily cache dataand restore the previous context."""

    # ==================================================================================
    def __init__(self, cache: LoaderCache) -> None:
        self._cache = cache
        self._token: Token[LoaderCache | None] | None = None

    # ==================================================================================
    def __enter__(self) -> LoaderCache:
        """Activate and return the session cache."""
        if self._token is not None:
            raise RuntimeError("A cache session cannot be re-entered.")

        self._token = LoaderCache._active_cache.set(self._cache)
        return self._cache

    # ==================================================================================
    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        """Clear the session and restore the previously active cache."""
        try:
            self._cache.clear()
        finally:
            if self._token is not None:
                LoaderCache._active_cache.reset(self._token)
                self._token = None

        return False
