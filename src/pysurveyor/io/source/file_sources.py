# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from dataclasses import dataclass

from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.source.abc.source import FileSource


# ======================================================================================
# CSV file source
# ======================================================================================
@dataclass
class CSVSource(FileSource[CSVConfig]):
    """Source containing a CSV file location and its :class:`CSVConfig`."""

    pass
