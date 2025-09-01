# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import numpy as np

from pysurv.validators._validators import validate_angle_unit
from pysurv.utils.utils import rho


def to_rad(angle: float, unit: str | None = None) -> float:
    """Convert angle from degrees or gradinas to radians."""
    unit = validate_angle_unit(unit)
    return angle / rho.get(unit, np.nan)


def from_rad(angle: float, unit: str | None = None) -> float:
    """Convert angle from radians to degrees or gradinas."""
    unit = validate_angle_unit(unit)
    return angle * rho.get(unit, np.nan)


def azimuth(x_first: float, y_first: float, x_second: float, y_second: float) -> float:
    """Calculate the azimuth in radians from coordinates."""
    # Convert inputs to arrays
    x_first = np.asarray(x_first)
    y_first = np.asarray(y_first)
    x_second = np.asarray(x_second)
    y_second = np.asarray(y_second)

    dx = x_second - x_first
    dy = y_second - y_first

    overlaps = (dx == 0) & (dy == 0)
    azimuths = np.full_like(overlaps, np.nan, dtype=float)
    azimuths[~overlaps] = np.mod(np.arctan2(dy[~overlaps], dx[~overlaps]), 2 * np.pi)

    if azimuths.size == 1:
        return float(azimuths)
    return azimuths
