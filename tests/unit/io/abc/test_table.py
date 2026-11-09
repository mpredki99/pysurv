# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests for the DataFrame-like :class:`pysurveyor.io.abc.Table` interface."""

from __future__ import annotations

import pandas as pd
from pandas.testing import assert_frame_equal, assert_index_equal, assert_series_equal
import pytest

from pysurveyor.io.abc.table import Table


# ======================================================================================
# Concrete Table implementation for exercising the abstract interface
# ======================================================================================
class TableProbe(Table[object]):  # type: ignore[misc]

    load = object()

    def copy(self, deep: bool = True):
        return TableProbe(self._data.copy(deep=deep))

    def _copy_kwargs(self):
        return {}


# ======================================================================================
# Pandas DataFrame fixture for testing the Table interface
# ======================================================================================
@pytest.fixture
def dataframe():
    return pd.DataFrame(
        {"name": ["A", "B", "C"], "value": [1.5, 2.5, 3.5]},
        index=pd.Index([10, 20, 30], name="record"),
    )


# ======================================================================================
# Table object fixture with a copy of the DataFrame
# ======================================================================================
@pytest.fixture
def table(dataframe):
    return TableProbe(dataframe.copy(deep=True))


# ======================================================================================
# Tests for the Table interface
# ======================================================================================
def test_item_access(dataframe, table):
    """Return the same columns and rows as the wrapped DataFrame."""
    assert_series_equal(table["value"], dataframe["value"])
    assert_frame_equal(table[["name", "value"]], dataframe[["name", "value"]])


# ======================================================================================
def test_item_assignment(dataframe, table):
    """Apply column assignments to the wrapped DataFrame in the same way."""
    dataframe["score"] = [10, 20, 30]
    table["score"] = [10, 20, 30]

    assert_frame_equal(table._data, dataframe)


# ======================================================================================
def test_properties(dataframe, table):
    """Expose the same index, columns, and dtypes as pandas."""
    assert_index_equal(table.index, dataframe.index)
    assert_index_equal(table.columns, dataframe.columns)
    assert_series_equal(table.dtypes, dataframe.dtypes)


# ======================================================================================
@pytest.mark.parametrize(
    ("accessor", "key"),
    [
        ("loc", (20, "value")),
        ("loc", (slice(10, 20), ["name", "value"])),
        ("iloc", (1, 1)),
        ("iloc", (slice(0, 2), [0, 1])),
        ("at", (30, "name")),
        ("iat", (2, 1)),
    ],
)
def test_indexers(dataframe, table, accessor, key):
    """Return the same values through each exposed pandas indexer."""
    table_result = getattr(table, accessor)[key]
    pandas_result = getattr(dataframe, accessor)[key]

    # ----------------------------------------------------------------------------------

    if isinstance(pandas_result, pd.DataFrame):
        assert_frame_equal(table_result, pandas_result)
    elif isinstance(pandas_result, pd.Series):
        assert_series_equal(table_result, pandas_result)
    else:
        assert table_result == pandas_result


# ======================================================================================
def test_equals(dataframe, table):
    """Compare equal and unequal data for Table and DataFrame."""
    equal_table = TableProbe(dataframe.copy(deep=True))

    different_dataframe = dataframe.copy(deep=True)
    different_dataframe.loc[10, "value"] = 99.0
    different_table = TableProbe(different_dataframe)

    # ----------------------------------------------------------------------------------

    assert table.equals(dataframe)
    assert table.equals(equal_table)
    assert not table.equals(different_dataframe)
    assert not table.equals(different_table)


# ======================================================================================
def test_equals_unsupported_operands(table):
    """Raise an error for values outside the supported interface."""
    with pytest.raises(TypeError):
        table.equals([1, 2, 3])  # type: ignore[arg-type]
