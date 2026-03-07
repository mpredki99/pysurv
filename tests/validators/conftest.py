# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def text_data():
    class TextData:
        empty = pd.Series([])
        valid = pd.Series(["a", "b", "c", "d", "e", "f"])
        convertible = pd.Series(["a", 1, pd.NA, "b", None, True])
        with_empty = pd.Series(["a", "b", pd.NA, "c", None, "d"])
        regex = pd.Series(
            [0, "A", 1, "# comment", "  # leading space", "	# leading tab"]
        )
        callable = pd.Series(["1", "A", "2", "B", "3", "C"])

    return TextData()


@pytest.fixture
def numeric_data():
    class NumericData:
        empty = pd.Series([])
        valid = pd.Series([-3, -2, -1, 0, 1, 2, 3])
        invalid = pd.Series(["A", "B", "C", "D", "E", "F", "G"])
        convertible = pd.Series(["-3", "-2", "-1", "0", "1", "2", "3"])
        with_empty = pd.Series([-3, 0, 3, pd.NA, np.nan, None, pd.NaT])
        regex = pd.Series([-3, -2, -1, "# 0", " #1", "	# 2", 3])
        callable = pd.Series([3, 2, 1, 0, -1, -2, -3])

    return NumericData()


@pytest.fixture
def equal_data():
    class NumericData:
        empty = pd.Series([])
        valid = pd.Series([0, 0, 0])
        invalid_type = pd.Series(["A", "B", "C"])
        invalid_value = pd.Series([0, 1, 2])
        with_empty = pd.Series([0, None, pd.NA])
        regex = pd.Series([-1, "# 0", " #1"])
        callable = pd.Series([-1, 0, 1])

    return NumericData()


@pytest.fixture
def comparison_data():
    class NumericData:
        empty = pd.Series([])
        valid = pd.Series([0, 1, 2])
        invalid_type = pd.Series(["A", "B", "C"])
        invalid_value = pd.Series([-1, 0, 1])
        with_empty = pd.Series([0, None, pd.NA])
        regex = pd.Series([0, "# 0", " #1"])
        callable = pd.Series([-2, 1, 3])

    return NumericData()


@pytest.fixture
def between_data():
    class BetweenData:
        empty = pd.Series([])
        valid = pd.Series([0, 100, 200, 300])
        invalid = pd.Series(["A", "B", "C", "D"])
        with_empty = pd.Series([0, 200, None, pd.NA])
        regex = pd.Series([0, 200, "# 300", " #400"])
        callable = pd.Series([0, 100, 111, 123])

    return BetweenData()
