# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from pysurveyor.backend.backend import DataFrame
from pysurveyor.io.accessors.loader_accessor import LoaderAccessor
from pysurveyor.io.accessors.parser_accessor import ParserAccessor
from pysurveyor.io.loader.schema_loader import SchemaLoader
from pysurveyor.io.loader.schema_parser.schema_parser import SchemaParser
from pysurveyor.schema.abc.strict_schema import StrictSchema
from pysurveyor.schema.model.models import models
from pysurveyor.schema.utils.validation_mode import ValidationMode


class PointsSchema(StrictSchema[SchemaLoader]):
    """
    Base class for points schema in the PySurveyor framework, providing
    a points schema representation with an associated validation model and column
    management.
    """

    __pysurveyor_abstract__ = True

    load = LoaderAccessor(SchemaLoader["PointsSchema"])
    parse = ParserAccessor(SchemaParser["PointsSchema"])

    # ----------------------------------------------------------------------------------
    _expected_index = models.points.index

    def __init__(
        self,
        model: DataFrame | None = None,
        validation_mode: ValidationMode | str = ValidationMode.LAZY,
    ) -> None:
        # ------------------------------------------------------------------------------
        # Initialize with the default model if no model is provided
        # ------------------------------------------------------------------------------
        super().__init__(
            models.points if model is None else model,
            validation_mode=validation_mode,
        )
