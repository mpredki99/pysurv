# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Unit tests for basic surveying calculations."""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from pysurveyor.basic.basic import (
    ang_subtract,
    azimuth,
    distance_3d,
    wrap_2pi,
    wrap_pi,
)

PI = np.pi


# ======================================================================================
# Unit tests for wrap_pi function
# ======================================================================================
class TestWrapPi:
    """Test wrapping angles to the ``[-pi, pi]`` interval."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("angle", "expected"),
        [
            (-2 * PI, 0.0),
            (-7 / 4 * PI, 1 / 4 * PI),
            (-3 / 2 * PI, 1 / 2 * PI),
            (-5 / 4 * PI, 3 / 4 * PI),
            (-PI, PI),
            (-3 / 4 * PI, -3 / 4 * PI),
            (-1 / 2 * PI, -1 / 2 * PI),
            (-1 / 4 * PI, -1 / 4 * PI),
            (0.0, 0.0),
            (1 / 4 * PI, 1 / 4 * PI),
            (1 / 2 * PI, 1 / 2 * PI),
            (3 / 4 * PI, 3 / 4 * PI),
            (PI, PI),
            (5 / 4 * PI, -3 / 4 * PI),
            (3 / 2 * PI, -1 / 2 * PI),
            (7 / 4 * PI, -1 / 4 * PI),
            (2 * PI, 0.0),
        ],
    )
    def test_scalar(self, angle, expected):
        """Test wrap_pi with scalar inputs."""
        result = wrap_pi(angle)

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(expected)

    # ==================================================================================
    def test_array(self):
        """Test wrap_pi with array input."""
        angles = np.array([0, PI / 2, PI, 3 * PI / 2, 2 * PI])
        result = wrap_pi(angles)

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert result.dtype == np.float64
        assert_allclose(result, [0.0, PI / 2, PI, -PI / 2, 0.0])

    # ==================================================================================
    def test_list_input(self):
        """Test wrap_pi with list input."""
        result = wrap_pi([0, PI, 2 * PI])

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert_allclose(result, [0.0, PI, 0.0])

    # ==================================================================================
    def test_zero_dimensional_array(self):
        """Test wrap_pi with zero-dimensional array input."""
        result = wrap_pi(np.array(PI))

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(PI)

    # ==================================================================================
    @pytest.mark.parametrize("periods", [1, 2, -1, -2])
    def test_periodicity(self, periods):
        """Test that wrap_pi is periodic with respect to multiples of 2 pi."""
        angles = np.arange(-3 * PI / 4, PI, PI / 4)
        result = wrap_pi(angles + periods * 2 * PI)

        # ------------------------------------------------------------------------------

        assert_allclose(result, angles)


# ======================================================================================
# Unit tests for wrap_2pi function
# ======================================================================================
class TestWrap2Pi:
    """Test wrapping angles to the ``[0, 2*pi)`` interval."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("angle", "expected"),
        [
            (-2 * PI, 0.0),
            (-7 / 4 * PI, 1 / 4 * PI),
            (-3 / 2 * PI, 1 / 2 * PI),
            (-5 / 4 * PI, 3 / 4 * PI),
            (-PI, PI),
            (-3 / 4 * PI, 5 / 4 * PI),
            (-PI / 2, 3 * PI / 2),
            (-PI / 4, 7 * PI / 4),
            (0.0, 0.0),
            (PI / 4, PI / 4),
            (PI / 2, PI / 2),
            (3 * PI / 4, 3 * PI / 4),
            (PI, PI),
            (5 * PI / 4, 5 * PI / 4),
            (3 * PI / 2, 3 * PI / 2),
            (7 * PI / 4, 7 * PI / 4),
            (2 * PI, 0.0),
        ],
    )
    def test_scalar(self, angle, expected):
        result = wrap_2pi(angle)

        assert isinstance(result, float)
        assert result == pytest.approx(expected)

    # ==================================================================================
    def test_array(self):
        angles = np.array([0, -PI / 2, -PI, -3 * PI / 2, -2 * PI])
        result = wrap_2pi(angles)

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert result.dtype == np.float64
        assert_allclose(result, [0.0, 3 * PI / 2, PI, PI / 2, 0.0])

    # ==================================================================================
    def test_list_input(self):
        result = wrap_2pi([0, PI, 2 * PI])

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert_allclose(result, [0.0, PI, 0.0])

    # ==================================================================================
    def test_zero_dimensional_array(self):
        result = wrap_2pi(np.array(PI))

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(PI)

    # ==================================================================================
    @pytest.mark.parametrize("periods", [1, 2, -1, -2])
    def test_periodicity(self, periods):
        angles = np.arange(0, 2 * PI, PI / 4)
        result = wrap_2pi(angles + periods * 2 * PI)

        # ------------------------------------------------------------------------------

        assert_allclose(result, angles)


# ======================================================================================
# Unit tests for ang_subtract function
# ======================================================================================
class TestAngSubtract:
    """Test circular subtraction of angles."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("a", "b", "expected"),
        [
            (1.0, 0.5, 0.5),
            (0.5, 1.0, -0.5),
            (PI, 0.0, PI),
            (0.0, PI, PI),
            (PI / 2, -PI / 2, PI),
            (0.0, 2 * PI, 0.0),
        ],
    )
    def test_scalar(self, a, b, expected):
        result = ang_subtract(a, b)

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(expected)

    # ==================================================================================
    def test_array(self):
        a = np.array([0.0, PI / 2, PI])
        b = np.array([PI / 2, 0.0, 0.0])

        result = ang_subtract(a, b)

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert_allclose(result, [-PI / 2, PI / 2, PI])

    # ==================================================================================
    def test_broadcasting(self):
        a = np.array([0.0])
        b = np.array([0.0, PI / 2])

        result = ang_subtract(a, b)

        # ------------------------------------------------------------------------------

        assert result.shape == (2,)
        assert_allclose(result, [0.0, -PI / 2])

    # ==================================================================================
    def test_list_inputs(self):
        result = ang_subtract([0, PI], [PI, 0])

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert_allclose(result, [PI, PI])

    # ==================================================================================
    def test_custom_angle_wrap(self):
        result = ang_subtract(0, PI / 4, angle_wrap=wrap_2pi)  # type: ignore

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(7 * PI / 4)

    # ==================================================================================
    def test_zero_dimensional_arrays(self):
        result = ang_subtract(np.array(1.0), np.array(0.5))

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(0.5)


# ======================================================================================
# Unit tests for azimuth function
# ======================================================================================
class TestAzimuth:
    """Test azimuth calculation between points."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("x1", "y1", "x2", "y2", "expected"),
        [
            (0, 0, 0, 1, 0.0),
            (0, 0, 1, 1, PI / 4),
            (0, 0, 1, 0, PI / 2),
            (0, 0, 1, -1, 3 * PI / 4),
            (0, 0, 0, -1, PI),
            (0, 0, -1, -1, 5 * PI / 4),
            (0, 0, -1, 0, 3 * PI / 2),
            (0, 0, -1, 1, 7 * PI / 4),
        ],
    )
    def test_scalar(self, x1, y1, x2, y2, expected):
        result = azimuth(x1, y1, x2, y2)

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(expected)

    # ==================================================================================
    def test_coincident_points(self):
        result = azimuth(1, 1, 1, 1)

        # ------------------------------------------------------------------------------

        assert np.isnan(result)

    # ==================================================================================
    def test_array(self):
        result = azimuth(
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 1, 0, -1],
            [1, 0, -1, 0],
        )

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert_allclose(result, [0.0, PI / 2, PI, 3 * PI / 2])

    # ==================================================================================
    def test_coincident_points_in_array(self):
        result = azimuth(
            [0, 1, 2],
            [0, 1, 2],
            [0, 1, 2],
            [0, 1, 2],
        )

        # ------------------------------------------------------------------------------

        assert_allclose(result, np.full(3, np.nan), equal_nan=True)

    # ==================================================================================
    def test_list_inputs(self):
        result = azimuth(
            [0, 0],
            [0, 0],
            [0, 1],
            [1, 0],
        )

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert_allclose(result, [0.0, PI / 2])

    # ==================================================================================
    def test_zero_dimensional_arrays(self):
        result = azimuth(
            np.array(0.0),
            np.array(0.0),
            np.array(1.0),
            np.array(0.0),
        )

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(PI / 2)

    # ==================================================================================
    def test_broadcasting(self):
        x1 = np.array([0.0])
        y1 = np.array([0.0])

        x2 = np.array([0.0, 1.0, -1.0])
        y2 = np.array([1.0, 1.0, 1.0])

        result = azimuth(x1, y1, x2, y2)

        # ------------------------------------------------------------------------------

        assert result.shape == (3,)
        assert_allclose(result, [0.0, PI / 4, 7 * PI / 4])


# ======================================================================================
# Unit tests for distance_3d function
# ======================================================================================
class TestDistance3D:
    """Test three-dimensional distance calculation."""

    # ==================================================================================
    def test_scalar(self):
        result = distance_3d(0, 0, 0, 1, 2, 2)

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(3.0)

    # ==================================================================================
    def test_coincident_points(self):
        result = distance_3d(1, 2, 3, 1, 2, 3)

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(0.0)

    # ==================================================================================
    def test_array(self):
        result = distance_3d(
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [1, 1, 0],
            [0, 2, 0],
            [0, 2, 3],
        )

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert result.dtype == np.float64
        assert_allclose(result, [1.0, 3.0, 3.0])

    # ==================================================================================
    def test_list_inputs(self):
        result = distance_3d(
            [0, 0],
            [0, 0],
            [0, 0],
            [3, 0],
            [4, 0],
            [0, 5],
        )

        # ------------------------------------------------------------------------------

        assert isinstance(result, np.ndarray)
        assert_allclose(result, [5.0, 5.0])

    # ==================================================================================
    def test_zero_dimensional_arrays(self):
        result = distance_3d(
            np.array(0.0),
            np.array(0.0),
            np.array(0.0),
            np.array(1.0),
            np.array(2.0),
            np.array(2.0),
        )

        # ------------------------------------------------------------------------------

        assert isinstance(result, float)
        assert result == pytest.approx(3.0)

    # ==================================================================================
    def test_broadcasting(self):
        result = distance_3d(
            0,
            0,
            0,
            np.array([1.0, 0.0, 0.0]),
            0,
            np.array([0.0, 2.0, 3.0]),
        )

        # ------------------------------------------------------------------------------

        assert result.shape == (3,)  # type: ignore
        assert_allclose(result, [1.0, 2.0, 3.0])
