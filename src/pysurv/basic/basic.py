# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import numpy as np
from numpy.typing import ArrayLike

from pysurv.typing.typing import AngleUnit
from pysurv.validators._validators import validate_angle_unit

from ._constants import RHO


def to_rad(angle: ArrayLike, unit: AngleUnit | None = None) -> ArrayLike:
    """
    Convert angle from degrees or gradinas to radians.

    Parameters
    ----------
    angle : ArrayLike
        Angle value(s) to convert.
    unit : AngleUnit, optional
        Unit of the input angle. Supported values are 'deg', 'grad', 'gon', or 'rad'.
        If None, a default unit is taken from pysurv.config object.

    Returns
    -------
    ArrayLike
        Angle value(s) in radians.

    Raises
    ------
    InvalidAngleUnitError
        If the unit parameter value is not supported.

    Examples
    --------
    >>> to_rad(180, 'deg')
    3.141592653589793

    >>> to_rad([200, 100], 'gon')
    array([3.14159265, 1.57079633])

    >>> to_rad(np.array([90, 180]), 'deg')
    array([1.57079633, 3.14159265])

    >>> to_rad(1.0, 'rad')
    1.0
    """
    unit = validate_angle_unit(unit)
    angles = np.asarray(angle)
    angles = angles / RHO.get(unit, np.nan)
    return angles.item() if angles.size == 1 else angles


def from_rad(angle: ArrayLike, unit: AngleUnit | None = None) -> ArrayLike:
    """
    Convert angle from radians to degrees or gradinas.

    Parameters
    ----------
    angle : ArrayLike
        Angle value(s) in radians to convert.
    unit : AngleUnit, optional
        Target unit for the output angle. Supported values are 'deg', 'grad', 'gon', or 'rad'.
        If None, a default unit is taken from pysurv.config object.

    Returns
    -------
    ArrayLike
        Angle value(s) in the specified unit.

    Raises
    ------
    InvalidAngleUnitError
        If the unit parameter value is not supported.

    Examples
    --------
    >>> from_rad(np.pi, 'deg')
    180.0

    >>> from_rad([np.pi, np.pi/2], 'gon')
    array([200., 100.])

    >>> from_rad(np.array([np.pi, np.pi/2]), 'gon')
    array([200., 100.])

    >>> from_rad(1.0, 'rad')
    1.0
    """
    unit = validate_angle_unit(unit)
    angles = np.asarray(angle)
    angles = angles * RHO.get(unit, np.nan)
    return angles.item() if angles.size == 1 else angles


def azimuth(
    x_first: ArrayLike, y_first: ArrayLike, x_second: ArrayLike, y_second: ArrayLike
) -> ArrayLike:
    """
    Calculate the azimuth in radians from coordinates.

    Azimuths are horizontal angles observed clockwise from any reference meridian.
    In plane surveying, azimuths are generally observed from north direction.

    Parameters
    ----------
    x_first : ArrayLike
        X coordinate(s) of the first point(s).
    y_first : ArrayLike
        Y coordinate(s) of the first point(s).
    x_second : ArrayLike
        X coordinate(s) of the second point(s).
    y_second : ArrayLike
        Y coordinate(s) of the second point(s).

    Returns
    -------
    ArrayLike
        Azimuth(s) in radians. If the input arrays are scalar, returns a scalar.
        If the points overlap, returns NaN for those positions.

    Examples
    --------
    >>> azimuth(0, 0, 1, 0)
    0.0

    >>> azimuth(0, 0, 0, 1)
    1.5707963267948966

    >>> azimuth(0, 0, -1, 0)
    3.141592653589793

    >>> azimuth(0, 0, 0, -1)
    4.71238898038469

    >>> azimuth([0, 0], [0, 0], [1, 0], [0, 1])
    array([0.,         1.57079633])

    >>> azimuth(1, 1, 1, 1)
    nan
    """
    # Convert inputs to arrays
    x_first, y_first, x_second, y_second = map(
        np.asarray, (x_first, y_first, x_second, y_second)
    )

    dx = x_second - x_first
    dy = y_second - y_first

    overlaps = (dx == 0) & (dy == 0)
    azimuths = np.full_like(overlaps, np.nan, dtype=float)
    azimuths[~overlaps] = np.mod(np.arctan2(dy[~overlaps], dx[~overlaps]), 2 * np.pi)

    return azimuths.item() if azimuths.size == 1 else azimuths
