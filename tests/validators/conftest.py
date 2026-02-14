import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def is_numeric_data():
    class NumericData:
        empty = pd.Series([])
        valid = pd.Series([-3, -2, -1, 0, 1, 2, 3])
        invalid = pd.Series(["A", "B", "C", "D", "E", "F", "G"])
        convertible = pd.Series(["-3", "-2", "-1", "0", "1", "2", "3"])
        with_empty = pd.Series([-3, 0, 3, pd.NA, np.nan, None, pd.NaT])
        regex_callable = pd.Series([-3, -2, -1, "# 0", " #1", "2", 3])

    return NumericData()


@pytest.fixture
def ge_data():
    class GeData:
        empty = pd.Series([])
        valid = pd.Series([0, 1, 2, 3, 4, 5, 6, 7])
        invalid = pd.Series([-3, -2, -1, 0, 1, 2, 3])
        with_empty = pd.Series([-3, 0, 3, pd.NA, np.nan, None, pd.NaT])

    return GeData()
