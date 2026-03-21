# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.


import pytest

from pysurv.schema import ControlsSchema

column_properties = [
    "coordinate_columns",
    "coordinate_sigma_columns",
]


@pytest.mark.parametrize("column_property", column_properties)
def test_columns_properties(column_property):
    schema = ControlsSchema()
    columns_set = getattr(schema, column_property)
    assert set(columns_set).issubset(schema.index)
