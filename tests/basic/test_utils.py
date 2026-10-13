"""Tests for `basic` module common utilities."""

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from pysurveyor.basic.utils import (
    as_array,
    as_scalar,
    is_scalar,
    normalize_numeric_type,
)


# ======================================================================================
# Test is_scalar
# ======================================================================================
@pytest.mark.parametrize(
    "value",
    [
        1,
        -1,
        1.5,
        0.0,
        np.int32(2),
        np.int64(2),
        np.float32(1.5),
        np.float64(1.5),
        np.array(1),
        np.array(1.5),
    ],
)
def test_scalar_values(value):
    assert is_scalar(value)


# ======================================================================================
@pytest.mark.parametrize(
    "value",
    [
        [1, 2],
        (1, 2),
        np.array([1]),
        np.array([[1]]),
        "1",
        None,
        {"value": 1},
    ],
)
def test_non_scalar_values(value):
    assert not is_scalar(value)


# ======================================================================================
# Test as_scalar
# ======================================================================================
@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (1, 1.0),
        (2.5, 2.5),
        (np.int64(3), 3.0),
        (np.float32(4.5), 4.5),
        (np.array(6), 6.0),
    ],
)
def test_scalar_conversion(value, expected):
    result = as_scalar(value)

    assert isinstance(result, float)
    assert result == pytest.approx(expected)


# ======================================================================================
def test_invalid_input_raises():
    with pytest.raises((TypeError, ValueError)):
        as_scalar("not a number")


# ======================================================================================
# Test as_array
# ======================================================================================
def test_list_to_float_array():
    result = as_array([1, 2, 3])

    assert isinstance(result, np.ndarray)
    assert result.dtype == np.float64
    assert_array_equal(result, np.array([1.0, 2.0, 3.0]))


# ======================================================================================
def test_integer_array_to_float_array():
    result = as_array(np.array([1, 2, 3], dtype=np.int32))

    assert result.dtype == np.float64
    assert_array_equal(result, [1.0, 2.0, 3.0])


# ======================================================================================
def test_scalar_input_produces_zero_dimensional_array():
    result = as_array(5)

    assert isinstance(result, np.ndarray)
    assert result.shape == ()
    assert result.dtype == np.float64
    assert result.item() == 5.0


# ======================================================================================
def test_empty_list():
    result = as_array([])

    assert result.shape == (0,)
    assert result.dtype == np.float64


# ======================================================================================
def test_invalid_input():
    with pytest.raises((TypeError, ValueError)):
        as_array("not a number")


# ======================================================================================
# Test normalize_numeric_type
# ======================================================================================
@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (1, 1.0),
        (2.5, 2.5),
        (np.int64(3), 3.0),
        (np.float32(4.5), 4.5),
        (np.array(6), 6.0),
    ],
)
def test_scalar_values_return_float(value, expected):
    result = normalize_numeric_type(value)

    assert isinstance(result, float)
    assert result == pytest.approx(expected)


# ======================================================================================
@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ([1, 2], [1.0, 2.0]),
        ((3, 4), [3.0, 4.0]),
        (np.array([1, 2], dtype=np.int32), [1.0, 2.0]),
        (np.array([[1, 2], [3, 4]]), [[1.0, 2.0], [3.0, 4.0]]),
    ],
)
def test_array_values_return_float_array(value, expected):
    result = normalize_numeric_type(value)

    assert isinstance(result, np.ndarray)
    assert result.dtype == np.float64
    assert_array_equal(result, expected)


# ======================================================================================
def test_invalid_input_raises():
    with pytest.raises((TypeError, ValueError)):
        normalize_numeric_type("not a number")
