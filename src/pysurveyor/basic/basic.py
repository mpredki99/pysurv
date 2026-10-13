"""Basic angular and coordinate utilities.

This module provides angle wrapping, angular subtraction, and azimuth calculations.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import overload

import numpy as np
from numpy.typing import ArrayLike

from .utils import FloatArray, as_array, as_scalar, is_scalar, normalize_numeric_type


# ======================================================================================
# Angle wrapping
# ======================================================================================
def _wrap(angle: object, lower: float, period: float) -> float | FloatArray:
    """Wrap angles into the half-open interval [lower, lower + period)."""
    return (normalize_numeric_type(angle) - lower) % period + lower


# ======================================================================================
@overload
def wrap_pi(angle: int | float | np.number) -> float: ...


@overload
def wrap_pi(angle: ArrayLike) -> FloatArray: ...


def wrap_pi(angle: object) -> float | FloatArray:
    """Wrap angles into the interval (-pi, pi]."""
    result = _wrap(angle, -np.pi, 2.0 * np.pi)

    # Convert the half-open interval [-pi, pi) to (-pi, pi].
    if is_scalar(result):
        value = as_scalar(result)
        return np.pi if value <= -np.pi else value

    return as_array(np.where(result <= -np.pi, np.pi, result))


# ======================================================================================
@overload
def wrap_2pi(angle: int | float | np.number) -> float: ...


@overload
def wrap_2pi(angle: ArrayLike) -> FloatArray: ...


def wrap_2pi(angle: object) -> float | FloatArray:
    """Wrap angles into the interval [0, 2*pi)."""
    return _wrap(angle, 0.0, 2.0 * np.pi)


# ======================================================================================
# Angular subtraction
# ======================================================================================
@overload
def ang_subtract(
    a: int | float | np.number,
    b: int | float | np.number,
    angle_wrap: Callable[[object], float | FloatArray] = wrap_pi,
) -> float: ...


@overload
def ang_subtract(
    a: ArrayLike,
    b: ArrayLike,
    angle_wrap: Callable[[object], float | FloatArray] = wrap_pi,
) -> FloatArray: ...


def ang_subtract(
    a: object,
    b: object,
    angle_wrap: Callable[[object], float | FloatArray] = wrap_pi,
) -> float | FloatArray:
    """Calculate the wrapped angular difference a - b.

    Parameters
    ----------
    a, b : object
        Angle or array of angles in radians.
    angle_wrap : callable, optional
        Function used to wrap the angular difference. Defaults to :func:`wrap_pi`.

    Returns
    -------
    float or FloatArray
        Wrapped angular difference(s).
    """
    difference = normalize_numeric_type(a) - normalize_numeric_type(b)
    return angle_wrap(difference)


# ======================================================================================
# Shared coordinate utilities
# ======================================================================================
def _all_scalar(*values: object) -> bool:
    """Return whether all inputs are scalar numeric values."""
    return all(is_scalar(value) for value in values)


# ======================================================================================
# Azimuth
# ======================================================================================
@overload
def azimuth(
    x_first: int | float | np.number,
    y_first: int | float | np.number,
    x_second: int | float | np.number,
    y_second: int | float | np.number,
) -> float: ...


@overload
def azimuth(
    x_first: ArrayLike,
    y_first: ArrayLike,
    x_second: ArrayLike,
    y_second: ArrayLike,
) -> FloatArray: ...


def azimuth(
    x_first: ArrayLike,
    y_first: ArrayLike,
    x_second: ArrayLike,
    y_second: ArrayLike,
) -> float | FloatArray:
    """Calculate the azimuth from the first point to the second.

    The azimuth is measured clockwise from the positive Y-axis and
    normalized to [0, 2*pi).

    Coincident points produce NaN.

    Parameters
    ----------
    x_first, y_first : ArrayLike
        Coordinates of the first point or points.
    x_second, y_second : ArrayLike
        Coordinates of the second point or points.

    Returns
    -------
    float or FloatArray
        Direction angle(s) in radians.

    Examples
    --------
    >>> azimuth(0, 0, 0, 1)
    0.0
    >>> azimuth(0, 0, 1, 0)
    1.5707963267948966
    >>> azimuth(0, 0, 0, -1)
    3.141592653589793
    >>> azimuth(0, 0, -1, 0)
    4.71238898038469
    >>> azimuth(1, 1, 1, 1)
    nan
    """
    scalar_input = _all_scalar(x_first, y_first, x_second, y_second)

    x_first, y_first, x_second, y_second = (
        as_array(x_first),
        as_array(y_first),
        as_array(x_second),
        as_array(y_second),
    )

    dx = x_second - x_first
    dy = y_second - y_first

    # Clockwise from +Y: east is pi/2, south is pi, west is 3*pi/2.
    result = np.mod(np.arctan2(dx, dy), 2.0 * np.pi)

    # A direction is undefined when the points coincide.
    result = np.where((dx == 0.0) & (dy == 0.0), np.nan, result)

    if scalar_input:
        return as_scalar(result)

    return as_array(result)


# ======================================================================================
# Euclidean distances
# ======================================================================================
@overload
def distance_3D(
    x_first: int | float | np.number,
    y_first: int | float | np.number,
    z_first: int | float | np.number,
    x_second: int | float | np.number,
    y_second: int | float | np.number,
    z_second: int | float | np.number,
) -> float: ...


@overload
def distance_3D(
    x_first: ArrayLike,
    y_first: ArrayLike,
    z_first: ArrayLike,
    x_second: ArrayLike,
    y_second: ArrayLike,
    z_second: ArrayLike,
) -> FloatArray: ...


def distance_3d(
    x_first: ArrayLike,
    y_first: ArrayLike,
    z_first: ArrayLike,
    x_second: ArrayLike,
    y_second: ArrayLike,
    z_second: ArrayLike,
) -> float | FloatArray:
    """Calculate the Euclidean distance between two 3D points.

    Coordinates support NumPy broadcasting. Scalar inputs return a Python
    float; array inputs return a float64 NumPy array.

    Parameters
    ----------
    x_first, y_first, z_first : ArrayLike
        Coordinates of the first point or points.
    x_second, y_second, z_second : ArrayLike
        Coordinates of the second point or points.

    Returns
    -------
    float or FloatArray
        Euclidean distance(s).

    Examples
    --------
    >>> distance_3D(0, 0, 0, 1, 2, 2)
    3.0
    """
    scalar_input = _all_scalar(x_first, y_first, z_first, x_second, y_second, z_second)

    x_first, y_first, z_first, x_second, y_second, z_second = (
        as_array(x_first),
        as_array(y_first),
        as_array(z_first),
        as_array(x_second),
        as_array(y_second),
        as_array(z_second),
    )

    dx = x_second - x_first
    dy = y_second - y_first
    dz = z_second - z_first

    result = np.hypot(np.hypot(dx, dy), dz)

    if scalar_input:
        return as_scalar(result)

    return as_array(result)
