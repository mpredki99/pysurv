# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
This module provides utility functions for basic surveying calculations.

It includes functions to calculate the azimuth angle and to convert angles
between radians, degrees, and gradians (gons).
"""
import numpy as np


def azimuth(first_point: dict, second_point: dict) -> float:
    """
    Calculate the azimuth angle between two points in a 2D plane, measured from the positive x-axis.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - first_point: (dict): A dictionary with 'x' and 'y' coordinates for the first point.
    - second_point: (dict): A dictionary with 'x' and 'y' coordinates for the second point.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - float: The azimuth angle in radians, ranging from 0 to 2π.

    --------------------------------------------------------------------------------------------------------------------
    Raises:

    - ValueError: If the first and second points overlap.
    """
    dx = second_point['x'] - first_point['x']
    dy = second_point['y'] - first_point['y']

    if dx == 0 and dy == 0:
        raise ValueError('Start and end points overlap. '
                         f'Azimuth for point is indefinite. dx = {dx}, dy = {dy}')

    az = np.arctan2(dy, dx)
    return az % (2 * np.pi)


def to_rad(angle: float, unit: str = 'grad') -> float:
    """
    Convert an angle to radians from either gradians (gons) or degrees.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - angle: (float): The angle to be converted.
    - unit: (str, optional): The unit of the angle ('grad', 'gon', or 'deg'). Defaults to 'grad'.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - float: The angle in radians.

    --------------------------------------------------------------------------------------------------------------------
    Raises:

    - ValueError: If the unit is not 'grad', 'gon', or 'deg'.
    """
    if unit in ['grad', 'gon']:
        return angle * np.pi / 200
    elif unit == 'deg':
        return angle * np.pi / 180
    else:
        raise ValueError('Invalid unit. Use "grad", "gon", or "deg".')


def from_rad(angle: float, unit: str = 'grad') -> float:
    """
    Convert an angle from radians to either gradians (gons) or degrees.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - angle: (float): The angle to be converted.
    - unit: (str, optional): The unit to convert the angle to ('grad', 'gon', or 'deg'). Defaults to 'grad'.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - float: The angle in the specified unit.

    --------------------------------------------------------------------------------------------------------------------
    Raises:

    - ValueError: If the unit is not 'grad', 'gon', or 'deg'.
    """
    if unit in ['grad', 'gon']:
        return angle * 200 / np.pi
    elif unit == 'deg':
        return angle * 180 / np.pi
    else:
        raise ValueError('Invalid unit. Use "grad", "gon", or "deg".')
