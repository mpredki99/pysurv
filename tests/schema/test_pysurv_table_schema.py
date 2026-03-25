# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest

from pysurv.schema import MeasurementsSchema


def test_assert_column_subtypes_in_schema_index():
    with pytest.raises(ValueError):
        MeasurementsSchema(
            angular_measurement_columns=["not_in_model", "not_in_model_2", "vz", "vh"]
        )


def test_assert_column_subtypes_unique_in_the_same_category():
    with pytest.raises(ValueError):
        MeasurementsSchema(angular_measurement_columns=["a", "a", "vz", "vh"])


def test_assert_column_subtypes_unique_in_the_other_category():
    with pytest.raises(ValueError):
        MeasurementsSchema(
            angular_measurement_columns=["a", "hz", "vz", "vh"],
            angular_measurement_sigma_columns=["a", "hz", "svz", "svh"],
        )
