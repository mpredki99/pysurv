"""Tests for the unit_conversion module."""

import numpy as np
import pytest

from pysurveyor.basic.unit_converter import convert


# ======================================================================================
# Test angle conversions
# ======================================================================================
@pytest.mark.parametrize(
    ("source", "target", "value", "expected"),
    [
        # ------------------------------------------------------------------------------
        ("rad", "rad", np.pi, np.pi),
        ("rad", "mrad", np.pi, np.pi * 1000),
        ("mrad", "rad", np.pi * 1000, np.pi),
        ("rad", "gon", np.pi, 200.0),
        ("gon", "rad", 200.0, np.pi),
        ("rad", "grad", np.pi, 200.0),
        ("grad", "rad", 200.0, np.pi),
        ("rad", "c", np.pi, 20000.0),
        ("c", "rad", 20000.0, np.pi),
        ("rad", "cc", np.pi, 2000000.0),
        ("cc", "rad", 2000000.0, np.pi),
        ("rad", "degree", np.pi, 180.0),
        ("degree", "rad", 180.0, np.pi),
        ("rad", "arcmin", np.pi, 10800.0),
        ("arcmin", "rad", 10800.0, np.pi),
        ("rad", "arcsec", np.pi, 648000.0),
        ("arcsec", "rad", 648000.0, np.pi),
        # ------------------------------------------------------------------------------
        ("degree", "mm_over_meter", 45.0, 1000.0),
        ("degree", "mm_over_meter", 1.0, 17.45506493),
        ("mm_over_meter", "gon", 1000.0, 50.0),
        ("mrad", "mm_over_meter", 0.01, 0.01),
        # ------------------------------------------------------------------------------
        ("degree", "gon", 180.0, 200.0),
        ("degree", "arcmin", 1.0, 60.0),
        ("degree", "arcsec", 1.0, 3600.0),
        ("arcmin", "arcsec", 1.0, 60.0),
        # ------------------------------------------------------------------------------
        ("gon", "degree", 200.0, 180.0),
        ("gon", "c", 1.0, 100.0),
        ("gon", "cc", 1.0, 10000.0),
        ("c", "cc", 1.0, 100.0),
        # ------------------------------------------------------------------------------
        ("rad", "gon", -np.pi, -200.0),
        ("degree", "rad", -180.0, -np.pi),
        ("degree", "arcmin", -1.0, -60.0),
        ("c", "cc", -1.0, -100.0),
    ],
)
def test_angle_scalar_conversions(source, target, value, expected):
    result = getattr(getattr(convert, source).to, target)(value)

    assert isinstance(result, float)
    assert result == pytest.approx(expected)


# ======================================================================================
@pytest.mark.parametrize(
    ("source", "target", "values", "expected"),
    [
        ("rad", "degree", [0.0, np.pi / 2, np.pi], [0.0, 90.0, 180.0]),
        ("gon", "degree", [0.0, 100.0, 200.0], [0.0, 90.0, 180.0]),
        ("mrad", "rad", [100.0, 500.0, 1000.0], [0.1, 0.5, 1.0]),
        ("degree", "arcmin", [1.0, 2.0, 3.0], [60.0, 120.0, 180.0]),
    ],
)
def test_angle_array_conversions(source, target, values, expected):
    result = getattr(getattr(convert, source).to, target)(values)

    assert isinstance(result, np.ndarray)
    np.testing.assert_allclose(result, expected)


# ======================================================================================
def test_angle_conversion_accepts_numpy_scalar():
    result = convert.degree.to.rad(np.float64(180.0))

    assert isinstance(result, float)
    assert result == pytest.approx(np.pi)


# ======================================================================================
def test_angle_conversion_accepts_zero_dimensional_array():
    value = np.array(180.0)

    result = convert.degree.to.rad(value)

    assert isinstance(result, float)
    assert result == pytest.approx(np.pi)


# ======================================================================================
# Degrees, minutes, seconds (DMS)
# ======================================================================================
def test_dms_to_degrees():
    result = convert.dms.to.degree([12.0, 30.0, 0.0])

    assert result == pytest.approx(12.5)


# ======================================================================================
def test_dms_to_radians():
    result = convert.dms.to.rad([180.0, 0.0, 0.0])

    assert result == pytest.approx(np.pi)


# ======================================================================================
def test_dms_to_degrees_with_negative_angle():
    result = convert.dms.to.degree([-12.0, 30.0, 0.0])

    assert result == pytest.approx(-12.5)


# ======================================================================================
def test_dms_to_degrees_array():
    values = np.array(
        [
            [0.0, 0.0, 0.0],
            [12.0, 30.0, 0.0],
            [45.0, 15.0, 0.0],
            [-10.0, 30.0, 0.0],
        ]
    )
    expected = np.array([0.0, 12.5, 45.25, -10.5])

    result = convert.dms.to.degree(values)

    np.testing.assert_allclose(result, expected)


# ======================================================================================
def test_degrees_to_dms():
    result = convert.degree.to.dms(12.5)

    np.testing.assert_allclose(result, [12.0, 30.0, 0.0], atol=1e-10)


# ======================================================================================
def test_degrees_to_dms_negative_angle():
    result = convert.degree.to.dms(-12.5)

    np.testing.assert_allclose(result, [-12.0, 30.0, 0.0], atol=1e-10)


# ======================================================================================
def test_degrees_to_dms_array():
    values = np.array([0.0, 12.5, 45.25, -90.5])
    expected = np.array(
        [
            [0.0, 0.0, 0.0],
            [12.0, 30.0, 0.0],
            [45.0, 15.0, 0.0],
            [-90.0, 30.0, 0.0],
        ]
    )

    result = convert.degree.to.dms(values)

    assert result.shape == (4, 3)
    np.testing.assert_allclose(result, expected, atol=1e-9)


# ======================================================================================
@pytest.mark.parametrize(
    "values",
    [
        12.0,
        [12.0],
        [12.0, 30.0],
        [[12.0, 30.0]],
        np.zeros((2, 2, 3)),
    ],
)
def test_dms_invalid_shapes(values):
    with pytest.raises(ValueError, match="shape \(3,\) or \(N, 3\)"):
        convert.dms.to.degree(values)


# ======================================================================================
@pytest.mark.parametrize(
    "values",
    [
        [12.0, -1.0, 0.0],
        [12.0, 0.0, -1.0],
        [[12.0, 30.0, -1.0]],
    ],
)
def test_dms_negative_minutes_or_seconds(values):
    with pytest.raises(ValueError, match="non-negative"):
        convert.dms.to.degree(values)


# ======================================================================================
# Length conversions
# ======================================================================================
@pytest.mark.parametrize(
    ("source", "target", "value", "expected"),
    [
        # ------------------------------------------------------------------------------
        ("kilometer", "meter", 1.0, 1000.0),
        ("meter", "kilometer", 1.0, 0.001),
        ("hectometer", "meter", 1.0, 100.0),
        ("meter", "hectometer", 1.0, 0.01),
        ("meter", "meter", 1.0, 1.0),
        ("decimeter", "meter", 1.0, 0.1),
        ("meter", "decimeter", 1.0, 10.0),
        ("centimeter", "meter", 1.0, 0.01),
        ("meter", "centimeter", 1.0, 100.0),
        ("millimeter", "meter", 1.0, 0.001),
        ("meter", "millimeter", 1.0, 1000.0),
        ("micrometer", "meter", 1.0, 0.000001),
        ("meter", "micrometer", 1.0, 1000000.0),
        # ------------------------------------------------------------------------------
        ("mile", "meter", 1.0, 1609.344),
        ("meter", "mile", 1.0, 0.0006213711922373339),
        ("yard", "meter", 1.0, 0.9144),
        ("meter", "yard", 1.0, 1.0936132983377078),
        ("feet", "meter", 1.0, 0.3048),
        ("meter", "feet", 1.0, 3.280839895013123),
        ("inch", "meter", 1.0, 0.0254),
        ("meter", "inch", 1.0, 39.37007874015748),
        ("thou", "meter", 1.0, 0.0000254),
        ("meter", "thou", 1.0, 39370.07874015748),
        # ------------------------------------------------------------------------------
        ("kilometer", "meter", -1.0, -1000.0),
        ("mile", "meter", -1.0, -1609.344),
        ("millimeter", "meter", -1.0, -0.001),
    ],
)
def test_length_scalar_conversions(source, target, value, expected):
    result = getattr(getattr(convert, source).to, target)(value)

    assert isinstance(result, float)
    assert result == pytest.approx(expected)


# ======================================================================================
@pytest.mark.parametrize(
    ("source", "target", "values", "expected"),
    [
        ("kilometer", "meter", [1.0, 2.0, 3.0], [1000.0, 2000.0, 3000.0]),
        ("meter", "centimeter", [0.0, 1.0, 2.5], [0.0, 100.0, 250.0]),
        ("mile", "feet", [1.0, 2.0], [5280.0, 10560.0]),
        ("inch", "millimeter", [1.0, 2.0], [25.4, 50.8]),
    ],
)
def test_length_array_conversions(source, target, values, expected):
    result = getattr(getattr(convert, source).to, target)(values)

    assert isinstance(result, np.ndarray)
    assert result.dtype == np.float64
    np.testing.assert_allclose(result, expected)


# ======================================================================================
def test_length_conversion_accepts_numpy_scalar():
    result = convert.kilometer.to.meter(np.float64(2.0))

    assert isinstance(result, float)
    assert result == pytest.approx(2000.0)


# ======================================================================================
def test_length_conversion_accepts_zero_dimensional_array():
    result = convert.meter.to.centimeter(np.array(2.5))

    assert isinstance(result, float)
    assert result == pytest.approx(250.0)


# ======================================================================================
# Invalid units
# ======================================================================================
@pytest.mark.parametrize(
    ("source", "target", "message"),
    [
        ("invalid", "rad", "Invalid unit: 'invalid'"),
        ("rad", "invalid", "Invalid angular unit: 'invalid'"),
        ("invalid", "meter", "Invalid unit: 'invalid'"),
        ("meter", "invalid", "Invalid length unit: 'invalid'"),
    ],
)
def test_conversion_invalid_units(source, target, message):
    with pytest.raises(AttributeError, match=message):
        getattr(getattr(convert, source).to, target)(0)


# ======================================================================================
@pytest.mark.parametrize(
    ("source", "target", "message"),
    [
        ("meter", "rad", "length to angular unit: 'rad'"),
        ("rad", "meter", "angle to length unit: 'meter'"),
    ],
)
def test_conversion_mixed_dimensions_units(source, target, message):
    with pytest.raises(AttributeError, match=message):
        getattr(getattr(convert, source).to, target)(0)


# =============================================================================
# Additional numerical edge cases
# =============================================================================
def test_empty_length_array():
    result = convert.meter.to.centimeter(np.array([]))

    assert result.shape == (0,)


def test_empty_angle_array():
    result = convert.degree.to.rad(np.array([]))

    assert result.shape == (0,)
