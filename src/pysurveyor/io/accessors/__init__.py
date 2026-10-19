# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Accessors for PySurveyor IO classes.

This package provides class descriptors that expose loaders and parsers
bound to their owning classes or objects.
"""

from .loader_accessor import LoaderAccessor
from .parser_accessor import ParserAccessor

__all__ = [
    "LoaderAccessor",
    "ParserAccessor",
]
