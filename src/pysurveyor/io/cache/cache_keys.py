# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any, Hashable, Mapping


# ======================================================================================
# File cache key
# ======================================================================================
@dataclass(frozen=True, slots=True)
class FileCacheKey:
    """Identify a file and the options used to read it.

    Attributes:
    ----------
    resolved_path : str
        Canonical absolute path to the file.
    reader_options : tuple[tuple[str, Hashable], ...]
        Sorted sequence of normalized key-value option pairs.
    modified_time_ns : int
        Last modification time of the file in nanoseconds.
    file_size : int
        Size of the file in bytes.
    """

    resolved_path: str
    reader_options: tuple[tuple[str, Hashable], ...]
    modified_time_ns: int
    file_size: int


# ======================================================================================
def create_file_cache_key(
    path: PathLike[str] | str, options: Mapping[str, Any]
) -> FileCacheKey:
    """Build a cache key from file metadata and reader options."""
    file_path = Path(path).expanduser().resolve()
    metadata = file_path.stat()

    # ----------------------------------------------------------------------------------
    reader_options = tuple(
        sorted((str(name), _make_hashable(value)) for name, value in options.items())
    )

    # ----------------------------------------------------------------------------------
    return FileCacheKey(
        resolved_path=str(file_path),
        reader_options=reader_options,
        modified_time_ns=metadata.st_mtime_ns,
        file_size=metadata.st_size,
    )


# ======================================================================================
# Internal utility to convert input values to hashable
# ======================================================================================
def _make_hashable(value: Any) -> Hashable:
    """Convert supported values into hashable cache-key components.

    Recursively transforms mappings, sequences, sets, and callables into canonical
    hashable types (such as tuples and frozensets).

    Parameters:
    ----------
    value : Any
        Value or data structure to convert into a hashable form.

    Returns:
    -------
    Hashable
        Hashable representation of the input value.

    Raises:
    ------
    TypeError
        If `value` or any nested element is not hashable and cannot be converted.
    """
    # Recursively convert mappings into hashable tuples --------------------------------
    if isinstance(value, Mapping):
        items = [
            (_make_hashable(key), _make_hashable(item)) for key, item in value.items()
        ]
        return ("mapping", tuple(sorted(items, key=repr)))

    # Recursively convert sequences into hashable tuples -------------------------------
    if isinstance(value, (list, tuple)):
        return ("sequence", tuple(_make_hashable(item) for item in value))

    # Recursively convert sets into hashable frozensets --------------------------------
    if isinstance(value, (set, frozenset)):
        return ("set", frozenset(_make_hashable(item) for item in value))

    # Convert callables into a hashable representation ---------------------------------
    if callable(value):
        module = getattr(value, "__module__", type(value).__module__)
        name = getattr(value, "__qualname__", repr(value))
        return ("callable", module, name)

    # Try to hash the value directly, raising an error if it is not hashable -----------
    try:
        hash(value)
    except TypeError as exc:
        raise TypeError(f"Unsupported cache-key value: {value!r}") from exc

    # ----------------------------------------------------------------------------------
    return value
