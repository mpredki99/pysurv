# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
The module contains functions for determining the coefficients of the observation equations of surveying measurements.
"""
from typing import Tuple
import numpy as np
from pysurv.modules.basic import azimuth


# Observation equation functions
def SD_obs_eq(meas_SD: float,
              idx_from: dict,
              idx_to: dict,
              coord_differences: dict,
              X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for slope distance (SD).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - meas_SD: (float): Measured slope distance.
    - idx_from: (dict): Indices for the station point:
        * x: (int): The index of the x-coordinate increment column of the station.
        * y: (int): The index of the y-coordinate increment column of the station.
        * z: (int): The index of the z-coordinate increment column of the station.
    - idx_to: (dict): Indices for the aim point:
        * x: (int): The index of the x-coordinate increment column of the aim point.
        * y: (int): The index of the y-coordinate increment column of the aim point.
        * z: (int): The index of the z-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dx: (float): The difference of the x-coordinates of the station and the aim point.
        * dy: (float): The difference of the y-coordinates of the station and the aim point.
        * dz: (float): The difference of the z-coordinates of the station and the aim point.
    - X_row: (np.ndarray): The X matrix row to be populated.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - tuple: A tuple of arrays:
        * (np.ndarray): updated X_row.
        * (np.ndarray): updated Y_row.
    """
    # Predefine necessary values
    distance = np.linalg.norm([coord_differences['dx'], coord_differences['dy'], coord_differences['dz']])
    # Calculate the coefficients of the matrix X
    X_row[idx_from['x']] = -np.divide(coord_differences['dx'], distance)
    X_row[idx_from['y']] = -np.divide(coord_differences['dy'], distance)
    X_row[idx_from['z']] = -np.divide(coord_differences['dz'], distance)
    X_row[idx_to['x']] = np.divide(coord_differences['dx'], distance)
    X_row[idx_to['y']] = np.divide(coord_differences['dy'], distance)
    X_row[idx_to['z']] = np.divide(coord_differences['dz'], distance)
    # Calculate the value of the free term
    Y_row = np.array(meas_SD - distance)

    return X_row, Y_row


def HD_obs_eq(meas_HD: float,
              idx_from: dict,
              idx_to: dict,
              coord_differences: dict,
              X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for horizontal distance (HD).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - meas_HD: (float): Measured horizontal distance.
    - idx_from: (dict): Indices for the station point:
        * x: (int): The index of the x-coordinate increment column of the station.
        * y: (int): The index of the y-coordinate increment column of the station.
    - idx_to: (dict): Indices for the aim point:
        * x: (int): The index of the x-coordinate increment column of the aim point.
        * y: (int): The index of the y-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dx: (float): The difference of the x-coordinates of the station and the aim point.
        * dy: (float): The difference of the y-coordinates of the station and the aim point.
    - X_row: (np.array): The X matrix row to be populated.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - tuple: A tuple of arrays:
        * (np.ndarray): updated X_row.
        * (np.ndarray): updated Y_row.
    """
    # Predefine necessary values
    distance = np.linalg.norm([coord_differences['dx'], coord_differences['dy']])
    # Calculate the coefficients of the matrix X
    X_row[idx_from['x']] = -np.divide(coord_differences['dx'], distance)
    X_row[idx_from['y']] = -np.divide(coord_differences['dy'], distance)
    X_row[idx_to['x']] = np.divide(coord_differences['dx'], distance)
    X_row[idx_to['y']] = np.divide(coord_differences['dy'], distance)
    # Calculate the value of the free term
    Y_row = np.array(meas_HD - distance)

    return X_row, Y_row


def dx_obs_eq(meas_dx: float,
              idx_from: dict,
              idx_to: dict,
              coord_differences: dict,
              X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for the x-coordinate component of the GNSS vector (dx).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - meas_dx: (float): Measured x-coordinate difference.
    - idx_from: (dict): Indices for the station point:
        * x: (int): The index of the x-coordinate increment column of the station.
    - idx_to: (dict): Indices for the aim point:
        * x: (int): The index of the x-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dx: (float): The difference of the x-coordinates of the station and the aim point.
    - X_row: (np.array): The X matrix row to be populated.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - tuple: A tuple of arrays:
        * (np.ndarray): updated X_row.
        * (np.ndarray): updated Y_row.
    """
    # Calculate the coefficients of the matrix X
    X_row[idx_from['x']] = -1
    X_row[idx_to['x']] = 1
    # Calculate the value of the free term
    Y_row = np.array(meas_dx - coord_differences['dx'])

    return X_row, Y_row


def dy_obs_eq(meas_dy: float,
              idx_from: dict,
              idx_to: dict,
              coord_differences: dict,
              X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for the y-coordinate component of the GNSS vector (dy).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - meas_dy: (float): Measured y-coordinate difference.
    - idx_from: (dict): Indices for the station point:
        * y: (int): The index of the y-coordinate increment column of the station.
    - idx_to: (dict): Indices for the aim point:
        * y: (int): The index of the y-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dy: (float): The difference of the y-coordinates of the station and the aim point.
    - X_row: (np.array): The X matrix row to be populated.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - tuple: A tuple of arrays:
        * (np.ndarray): updated X_row.
        * (np.ndarray): updated Y_row.
    """
    # Calculate the coefficients of the matrix X
    X_row[idx_from['y']] = -1
    X_row[idx_to['y']] = 1
    # Calculate the value of the free term
    Y_row = np.array(meas_dy - coord_differences['dy'])

    return X_row, Y_row


def dz_obs_eq(meas_dz: float,
              idx_from: dict,
              idx_to: dict,
              coord_differences: dict,
              X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for the z-coordinate component of the GNSS vector (dz).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - meas_dz: (float): Measured z-coordinate difference.
    - idx_from: (dict): Indices for the station point:
        * z: (int): The index of the z-coordinate increment column of the station.
    - idx_to: (dict): Indices for the aim point:
        * z: (int): The index of the z-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dz: (float): The difference of the z-coordinates of the station and the aim point.
    - X_row: (np.array): The X matrix row to be populated.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - tuple: A tuple of arrays:
        * (np.ndarray): updated X_row.
        * (np.ndarray): updated Y_row.
    """
    # Calculate the coefficients of the matrix X
    X_row[idx_from['z']] = -1
    X_row[idx_to['z']] = 1
    # Calculate the value of the free term
    Y_row = np.array(meas_dz - coord_differences['dz'])

    return X_row, Y_row


def A_obs_eq(meas_A: float,
             idx_from: dict,
             idx_to: dict,
             coord_differences: dict,
             X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for azimuth angle (A).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - meas_A: (float): Measured azimuth angle.
    - idx_from: (dict): Indices for the station point:
        * x: (int): The index of the x-coordinate increment column of the station.
        * y: (int): The index of the y-coordinate increment column of the station.
    - idx_to: (dict): Indices for the aim point:
        * x: (int): The index of the x-coordinate increment column of the aim point.
        * y: (int): The index of the y-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dx: (float): The difference of the x-coordinates of the station and the aim point.
        * dy: (float): The difference of the y-coordinates of the station and the aim point.
    - X_row: (np.array): The X matrix row to be populated.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - tuple: A tuple of arrays:
        * (np.ndarray): updated X_row.
        * (np.ndarray): updated Y_row.
    """
    # Predefine necessary values
    distance_sq = np.linalg.norm([coord_differences['dx'], coord_differences['dy']]) ** 2
    azimuth_0 = azimuth({'x': 0, 'y': 0}, {'x': coord_differences['dx'], 'y': coord_differences['dy']})
    # Calculate the coefficients of the matrix X
    X_row[idx_from['x']] = np.divide(coord_differences['dy'], distance_sq)
    X_row[idx_from['y']] = -np.divide(coord_differences['dx'], distance_sq)
    X_row[idx_to['x']] = -np.divide(coord_differences['dy'], distance_sq)
    X_row[idx_to['y']] = np.divide(coord_differences['dx'], distance_sq)
    # Calculate the value of the free term
    Y_row = np.array(meas_A - azimuth_0)

    return X_row, Y_row


def HZ_obs_eq(meas_HZ: float,
              idx_from: dict,
              idx_to: dict,
              coord_differences: dict,
              X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for horizontal direction (HZ).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - meas_HZ: (float): Measured horizontal direction.
    - idx_from: (dict): Indices for the station point:
        * x: (int): The index of the x-coordinate increment column of the station.
        * y: (int): The index of the y-coordinate increment column of the station.
        * o: (int): The index of the orientation constant column of the station.
    - idx_to: (dict): Indices for the aim point:
        * x: (int): The index of the x-coordinate increment column of the aim point.
        * y: (int): The index of the y-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dx: (float): The difference of the x-coordinates of the station and the aim point.
        * dy: (float): The difference of the y-coordinates of the station and the aim point.
        * o: (float): The orientation constant value of the station.
    - X_row: (np.array): The X matrix row to be populated.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - tuple: A tuple of arrays:
        * (np.ndarray): updated X_row.
        * (np.ndarray): updated Y_row.
    """
    # Predefine necessary values
    distance_sq = np.linalg.norm([coord_differences['dx'], coord_differences['dy']])**2
    azimuth_0 = azimuth({'x': 0, 'y': 0}, {'x': coord_differences['dx'], 'y': coord_differences['dy']})
    # Calculate the coefficients of the matrix X
    X_row[idx_from['x']] = np.divide(coord_differences['dy'], distance_sq)
    X_row[idx_from['y']] = -np.divide(coord_differences['dx'], distance_sq)
    X_row[idx_from['o']] = -1
    X_row[idx_to['x']] = - np.divide(coord_differences['dy'], distance_sq)
    X_row[idx_to['y']] = np.divide(coord_differences['dx'], distance_sq)
    # Calculate the value of the free term
    Y_row = np.array(meas_HZ - azimuth_0 + coord_differences['o'])
    # Reduce the value of the free term if necessary
    if Y_row >= (2 * np.pi - 0.02):
        Y_row -= np.pi * 2
    elif Y_row <= -(2 * np.pi - 0.02):
        Y_row += np.pi * 2

    return X_row, Y_row


def VZ_obs_eq(meas_VZ: float,
              idx_from: dict,
              idx_to: dict,
              coord_differences: dict,
              X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for vertical zenith angle (VZ).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - meas_VZ: (float): Measured vertical zenith angle.
    - idx_from: (dict): Indices for the station point:
        * x: (int): The index of the x-coordinate increment column of the station.
        * y: (int): The index of the y-coordinate increment column of the station.
        * z: (int): The index of the z-coordinate increment column of the station.
    - idx_to: (dict): Indices for the aim point:
        * x: (int): The index of the x-coordinate increment column of the aim point.
        * y: (int): The index of the y-coordinate increment column of the aim point.
        * z: (int): The index of the z-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dx: (float): The difference of the x-coordinates of the station and the aim point.
        * dy: (float): The difference of the y-coordinates of the station and the aim point.
        * dz: (float): The difference of the z-coordinates of the station and the aim point.

    - X_row: (np.array): The X matrix row to be populated.

    ------------
    Returns:

    - tuple: A tuple of arrays:
        - (np.ndarray): updated X_row.
        - (np.ndarray): updated Y_row.
    """
    # Predefine necessary values
    distance_HD = np.linalg.norm([coord_differences['dx'], coord_differences['dy']])
    distance_SD_sq = np.linalg.norm([coord_differences['dx'], coord_differences['dy'], coord_differences['dz']])**2
    # Calculate the coefficients of the matrix X
    X_row[idx_from['x']] = -np.divide(coord_differences['dx'] * coord_differences['dz'], distance_HD * distance_SD_sq)
    X_row[idx_from['y']] = -np.divide(coord_differences['dy'] * coord_differences['dz'], distance_HD * distance_SD_sq)
    X_row[idx_from['z']] = np.divide(distance_HD, distance_SD_sq)
    X_row[idx_to['x']] = np.divide(coord_differences['dx'] * coord_differences['dz'], distance_HD * distance_SD_sq)
    X_row[idx_to['y']] = np.divide(coord_differences['dy'] * coord_differences['dz'], distance_HD * distance_SD_sq)
    X_row[idx_to['z']] = -np.divide(distance_HD, distance_SD_sq)
    # Calculate the value of the free term
    Y_row = np.array(
        meas_VZ - np.arctan2(
            distance_HD, coord_differences['dz']
        )
    )
    return X_row, Y_row


def VH_obs_eq(meas_VH: float,
              idx_from: dict,
              idx_to: dict,
              coord_differences: dict,
              X_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the observation equation for vertical horizontal angle (VH).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - VH_obs_eq: (float): Measured vertical horizontal angle.
    - idx_from: (dict): Indices for the station point:
        * x: (int): The index of the x-coordinate increment column of the station.
        * y: (int): The index of the y-coordinate increment column of the station.
        * z: (int): The index of the z-coordinate increment column of the station.
    - idx_to: (dict): Indices for the aim point:
        * x: (int): The index of the x-coordinate increment column of the aim point.
        * y: (int): The index of the y-coordinate increment column of the aim point.
        * z: (int): The index of the z-coordinate increment column of the aim point.
    - coord_differences: (dict): Coordinate differences between the points:
        * dx: (float): The difference of the x-coordinates of the station and the aim point.
        * dy: (float): The difference of the y-coordinates of the station and the aim point.
        * dz: (float): The difference of the z-coordinates of the station and the aim point.
    - X_row: (np.array): The X matrix row to be populated.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - tuple: A tuple of arrays:
        * (np.ndarray): updated X_row.
        * (np.ndarray): updated Y_row.
    """
    # Predefine necessary values
    distance_HD = np.linalg.norm([coord_differences['dx'], coord_differences['dy']])
    distance_SD_sq = np.linalg.norm([coord_differences['dx'], coord_differences['dy'], coord_differences['dz']])**2
    # Calculate the coefficients of the matrix X
    X_row[idx_from['x']] = np.divide(coord_differences['dx'] * coord_differences['dz'], distance_HD * distance_SD_sq)
    X_row[idx_from['y']] = np.divide(coord_differences['dy'] * coord_differences['dz'], distance_HD * distance_SD_sq)
    X_row[idx_from['z']] = -np.divide(distance_HD, distance_SD_sq)
    X_row[idx_to['x']] = -np.divide(coord_differences['dx'] * coord_differences['dz'], distance_HD * distance_SD_sq)
    X_row[idx_to['y']] = -np.divide(coord_differences['dy'] * coord_differences['dz'], distance_HD * distance_SD_sq)
    X_row[idx_to['z']] = np.divide(distance_HD, distance_SD_sq)
    # Calculate the value of the free term
    Y_row = np.array(
        meas_VH - np.arctan2(coord_differences['dz'], distance_HD)
    )
    return X_row, Y_row


# Dictionary mapping observation types to their respective functions
observation_functions = {
    'SD': SD_obs_eq,
    'HD': HD_obs_eq,
    'VD': dz_obs_eq,
    'dx': dx_obs_eq,
    'dy': dy_obs_eq,
    'dz': dz_obs_eq,
    'A': A_obs_eq,
    'HZ': HZ_obs_eq,
    'VZ': VZ_obs_eq,
    'VH': VH_obs_eq,
}
