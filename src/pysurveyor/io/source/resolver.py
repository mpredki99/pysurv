# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Any

from pysurveyor.io.config.abc.config import IOConfig
from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.source.abc.source import Source
from pysurveyor.io.source.file_sources import CSVSource, FileSource


# ======================================================================================
# Public source resolver
# ======================================================================================
def resolve_source(
    input: Any,
    config: IOConfig | None = None,
) -> Source:
    """Automatically resolve an input into a concrete source."""
    if isinstance(input, (str, PathLike)):
        return resolve_file_source(Path(input), config)

    # ----------------------------------------------------------------------------------
    raise TypeError(f"Could not resolve source from input: {input!r}")


# ======================================================================================
# Source resolving helper function
# ======================================================================================
def resolve_file_source(path: Path, config: IOConfig | None = None) -> FileSource:
    """Resolve an existing file path into a format-specific source.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not identify an existing file.
    TypeError
        If the file suffix or configuration type is unsupported.
    """
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"File source does not exist: {path}")

    # ----------------------------------------------------------------------------------
    suffix = path.suffix.lower()

    match suffix:
        # ------------------------------------------------------------------------------
        case ".csv":
            if config is not None and not isinstance(config, CSVConfig):
                raise TypeError(
                    f"CSV sources require CSVConfig, got {type(config).__name__}"
                )
            return CSVSource(path=path, config=config or CSVConfig())
        # ------------------------------------------------------------------------------
        case _:
            raise TypeError(f"Unsupported file suffix: {suffix!r} for file: {path!r}")
