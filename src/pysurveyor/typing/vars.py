# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Shared type-variable definitions for generic PySurveyor classes.

All bounds use forward-reference strings to avoid import cycles at runtime.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

if TYPE_CHECKING:
    from pysurveyor.abc.pysurveyor_abc import PySurveyorABC
    from pysurveyor.data.abc.datatable import DataTable
    from pysurveyor.io.config.abc.config import IOConfig
    from pysurveyor.io.loader.abc.pysurveyor_loader import PySurveyorLoader
    from pysurveyor.io.loader.schema_parser.schema_parser import SchemaParser
    from pysurveyor.schema.abc.schema import Schema


# ======================================================================================
LoaderType = TypeVar("LoaderType")
SchemaParserType = TypeVar("SchemaParserType")

# ======================================================================================
DomainClassType = TypeVar("DomainClassType", bound="PySurveyorABC")

# ======================================================================================
ConfigType = TypeVar("ConfigType", bound="IOConfig")
