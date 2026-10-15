# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""PySurveyor data source classes and resolution helpers."""

from . import abc
from .file_sources import CSVSource
from .resolver import resolve_source

__all__ = [
    "abc",
    "CSVSource",
    "resolve_source",
]
