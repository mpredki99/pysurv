# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.


import pytest

from pysurv.schema import MeasurementsSchema

column_properties = [
    "target_columns",
    "linear_measurement_columns",
    "linear_measurement_sigma_columns",
    "angular_measurement_columns",
    "angular_measurement_sigma_columns",
    "measurement_columns",
    "measurement_sigma_columns",
    "linear_columns",
    "angular_columns",
]


@pytest.mark.parametrize("column_property", column_properties)
def test_columns_properties(column_property):
    schema = MeasurementsSchema()
    columns_set = getattr(schema, column_property)
    assert set(columns_set).issubset(schema.index)
