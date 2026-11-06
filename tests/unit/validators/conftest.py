# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Shared sample Series for validator tests."""

from types import SimpleNamespace

import numpy as np
from pandas import NA, Series
import pytest


# ======================================================================================
def _empty() -> Series:
    """Return an empty Series."""
    return Series([], dtype=object)


# ======================================================================================
# Sample data for conversion validators
# ======================================================================================
@pytest.fixture
def numeric_data() -> SimpleNamespace:
    """Provide numeric samples for ``IsNumeric`` tests."""
    return SimpleNamespace(
        empty=_empty(),
        valid=Series([-3, -2, -1, 0, 1, 2, 3]),
        invalid=Series(["A", "B", "C", "D", "E", "F", "G"]),
        convertible=Series(["-3", "-2", "-1", "0", "1", "2", "3"]),
        with_empty=Series([-3, 0, 3, NA, np.nan, None]),
        literal=Series(["A", False, -1, 0, 1, 2, 3]),
        regex=Series([-3, -2, -1, "# 0", " #1", "\t# 2", 3]),
        callable=Series([3, 2, 1, 0, -1, -2, -3]),
    )


# ======================================================================================
@pytest.fixture
def text_data() -> SimpleNamespace:
    """Provide text samples for ``IsText`` tests."""
    return SimpleNamespace(
        empty=_empty(),
        valid=Series(["a", "b", "c", "d", "e", "f"]),
        convertible=Series(["a", 1, NA, "b", None, True]),
        with_empty=Series(["a", "b", NA, "c", None, "d"]),
        literal=Series([1, True, "c", "d", "e", "f"]),
        regex=Series([0, "A", 1, "# comment", "  # leading space", "\t# leading tab"]),
        callable=Series(["1", "A", "2", "B", "3", "C"]),
    )


# ======================================================================================
@pytest.fixture
def boolean_data() -> SimpleNamespace:
    """Provide boolean samples for ``IsBoolean`` tests."""
    return SimpleNamespace(
        empty=_empty(),
        valid=Series([True, False, True, False, True, False]),
        invalid=Series(["A", "B", "C", "D", 1.01, 3.5, 0.99]),
        convertible=Series([1, 0, True, False, 1.0, 0.0]),
        with_empty=Series([True, None, NA, False, True, False]),
        literal=Series(["A", 1, True, False, True, False]),
        regex=Series([True, False, "# comment", "  # leading space", "\t# Tab", True]),
        callable=Series([10, 100, True, False, True, False]),
    )


# ======================================================================================
# Sample data for comparison validators
# ======================================================================================
@pytest.fixture
def equal_data() -> SimpleNamespace:
    """Provide equality-comparison samples."""
    return SimpleNamespace(
        empty=_empty(),
        valid=Series([0, 0, 0]),
        invalid_type=Series(["A", "B", "C"]),
        invalid_value=Series([0, 1, 2]),
        with_empty=Series([0, None, NA]),
        literal=Series([0, 1, 2]),
        regex=Series([-1, "# 0", " #1"]),
        callable=Series([-1, 0, 1]),
    )


# ======================================================================================
@pytest.fixture
def comparison_data() -> SimpleNamespace:
    """Provide inequality-comparison samples."""
    return SimpleNamespace(
        empty=_empty(),
        valid=Series([0, 1, 2]),
        invalid_type=Series(["A", "B", "C"]),
        invalid_value=Series([-1, 0, 1]),
        with_empty=Series([0, None, NA]),
        literal=Series([-1, 0, 1]),
        regex=Series([0, "# 0", " #1"]),
        callable=Series([-2, 1, 3]),
    )


# ======================================================================================
# Sample data for Between validator
# ======================================================================================
@pytest.fixture
def between_data() -> SimpleNamespace:
    """Provide range-comparison samples."""
    return SimpleNamespace(
        empty=_empty(),
        valid=Series([0, 100, 200, 300]),
        invalid=Series(["A", "B", "C", "D"]),
        with_empty=Series([0, 200, None, NA]),
        literal=Series([-100, 0, 100, 200]),
        regex=Series([0, 200, "# 300", " #400"]),
        callable=Series([0, 100, 111, 123]),
    )
