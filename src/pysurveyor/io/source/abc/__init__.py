# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Abstract base classes for PySurveyor data sources."""

from .source import FileSource, Source

__all__ = ["Source", "FileSource"]
