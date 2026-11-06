# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests of schema loader."""

import pytest

from pysurveyor.io.loader import SchemaLoader
from pysurveyor.schema import (
    ControlPointsSchema,
    InstrumentsSchema,
    MeasurementsSchema,
    PointsSchema,
    StationsSchema,
    TargetsSchema,
)


# ======================================================================================
@pytest.mark.parametrize(
    ("constructor"),
    [
        ControlPointsSchema,
        InstrumentsSchema,
        MeasurementsSchema,
        PointsSchema,
        StationsSchema,
        TargetsSchema,
    ],
)
def test_bound_constructor(constructor):
    """Test that the loader is bound to the correct constructor."""

    loader = SchemaLoader(constructor)

    # ----------------------------------------------------------------------------------

    assert loader._constructor is constructor
