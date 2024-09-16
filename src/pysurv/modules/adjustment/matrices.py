# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
A module for creating and managing matrices used in the adjustment of surveying control networks.

It supports the creation of observation equation matrices, weight matrices, and constraint matrices
for free adjustment, as well as matrices for adjustment including the standard deviations of control points.
"""
from typing import Tuple
import numpy as np
import pandas as pd
from pysurv.modules import Controls
from pysurv.modules import Measurements
from pysurv.modules.basic import azimuth
from pysurv import customizations
from .observation_equations import observation_functions


def equations_system(controls: Controls, measurements: Measurements,
                     methods: dict = customizations.methods['default'],
                     default_sigma: dict = customizations.default_measurement_sigma['default']) -> dict:
    """
    Creates a system of observation equations for adjustment of surveying control network.

    System of equations is represented as the matrices.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - controls: (Controls): The control points dataset.
    - measurements: (Measurements): The measurement dataset.
    - methods: (dict, optional): Dictionary with methods to use to set observation and control point weights,
      and also with tuning constants. Defaults as 'default' in customisations.
    - default_sigma: (dict, optional): Default standard deviations for measurements and point locations.
      Defaults as 'default' in customisations.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - dict: A dictionary containing matrices:
        - X: (np.ndarray): The coefficient matrix.
        - Y: (np.ndarray): The vector of observed values.
        - W: (np.ndarray, optional): The observations weight matrix.
        - R: (np.ndarray, optional): The matrix of network inner constraints.
        - sX: (np.ndarray, optional): The coefficient matrix for control points increments.
        - sW: (np.ndarray, optional): The weight matrix for the control points.
    """
    # Calculate approximate values of orientations if necessary
    if 'HZ' in measurements.types:
        controls['o'] = approx_orientation(controls, measurements['HZ'])

    # Create matrix X, Y and W if necessary
    matrices = XYW_matrices(
        controls,
        measurements,
        calculate_W=methods['observations'] != 'ordinary',
        default_sigma=default_sigma,
    )
    # Determine if matrix of inner constraint is needed
    free_method = methods.get('free')
    if free_method == 'ordinary':
        return matrices
    # Add matrix R
    elif free_method:
        matrices.update(
            apply_inner_constraints(
                controls,
                measurements.types
            )
        )
    # Calculate matrix sW and sX if needed
    matrices.update(
        calculate_sW(
            controls,
            calculate_sX=free_method is None,
            default_pt_sigma=default_sigma['sP']
        )
    )
    return matrices


def XYW_matrices(controls: Controls, measurements: Measurements, calculate_W: bool = True,
                 default_sigma: dict = customizations.default_measurement_sigma['default']) -> dict:
    """
    Creates a system of observation equations for surveying measurements.

    This function iterates through each pair of points specified in the measurements' dataset,
    calculates the necessary observation equations using the appropriate observation function,
    and builds the X (coefficient matrix), Y (observation vector), and optionally W (diagonal weight matrix) matrices.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - controls: (Controls): Coordinates of control points.
    - measurements: (Measurements): The measurement dataset.
    - calculate_W: (bool, optional): If True, the function calculates the weight matrix (W) for the observations.
      Defaults to True.
    - default_sigma: (dict, optional): Default standard deviations for measurements and point locations.
      Defaults as 'default' in customisations.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - dict: A dictionary containing matrices:
        - X: (np.ndarray): The coefficient matrix.
        - Y: (np.ndarray): The vector of observed values.
        - W: (np.ndarray, optional): The observations weight matrix.
    """
    # Cache variables
    controls_coords = controls[controls.labels]
    coord_indices = controls.coord_index
    n_measurements = measurements.n_measurements
    n_coords = controls.n_coords
    measuremnts_types = measurements.types
    # Initialize matrices
    X = np.zeros((n_measurements, n_coords))
    Y = np.zeros((n_measurements, 1))
    W = np.zeros((n_measurements, n_measurements)) if calculate_W else None
    # Measurements counter
    i_meas = 0
    for (meas_from, meas_to), row in measurements.iterrows():
        # Determine column indices in matrix X for control points labels
        idx_from, idx_to = coord_indices[meas_from], coord_indices[meas_to]
        # Calculate coordinate differences between the control points 'FROM' and 'TO'
        coord_differences = _get_coords_differences(controls_coords.loc[meas_from], controls_coords.loc[meas_to])

        for measurement_type in measuremnts_types:
            measurement_value = row[measurement_type]
            # Skip the measurement if it is empty
            if np.isnan(measurement_value):
                continue
            # Calculate observation equation coefficients for the current measurement
            X[i_meas, :], Y[i_meas] = observation_functions[measurement_type](
                measurement_value,
                idx_from,
                idx_to,
                coord_differences,
                X[i_meas, :]
            )
            # Calculate the weight for the current measurement if necessary
            if W is not None:
                W[i_meas, i_meas] = set_obs_weight(
                    row,
                    measurement_type,
                    measurement_value,
                    default_sigma
                )
            # Increase measurements counter
            i_meas += 1
    # Save results in dictionary
    matrices = {'X': X, 'Y': Y}
    if W is not None:
        matrices['W'] = W

    return matrices


def apply_inner_constraints(controls: Controls, measurement_types: pd.Index) -> dict:
    """
    Applies inner constraints for free network adjustment.

    Function creates the R matrix (matrix of inner constraints).

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - controls: (Controls): Coordinates of control points.
    - measurement_types: (pd.Index): A list of measurement types in dataset.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - dict: A dictionary with 'R' matrix and list of applied constraints:
        - R: (np.ndarray): Matrix of inner constraints.
        - constraints: (list): List of applied constraints.
    """
    # Initialize matrix and list of applied constraints
    R = list()
    constraints_list = list()
    # Get indices of columns with appropriate coordinate label
    coord_idx = controls.coord_index
    coord_indices = {label:
        [coord_idx[point][label] for point in coord_idx if label in coord_idx[point]]
        for label in controls.labels
    }
    # Create dict of observations conditions
    obs_condition = {
        '1D': not ('x' in controls.labels and 'y' in controls.labels) and 'z' in controls.labels,
        '3D': all(coord in controls.labels for coord in ['x', 'y', 'z']),
        'SD': 'SD' in measurement_types or
              ('HD' in measurement_types and 'VD' in measurement_types) or
              ('HD' in measurement_types and 'VZ' in measurement_types) or
              ('VD' in measurement_types and 'VZ' in measurement_types) or
              ('HD' in measurement_types and 'VH' in measurement_types) or
              ('VD' in measurement_types and 'VH' in measurement_types),
        'HD': 'HD' in measurement_types,
        'VD': any(measurement in measurement_types for measurement in ['VD', 'dz']),
        'vector': all(measurement in measurement_types for measurement in ['dx', 'dy', 'dz']),
        'HZ': 'HZ' in measurement_types,
        'VZ': 'VZ' in measurement_types,
        'A': 'A' in measurement_types,
    }
    # Create dict of constraints conditions
    constraint_condition = {
        'rotate': not obs_condition['1D'] and not (obs_condition['vector'] or obs_condition['A']),
        '1D_scale': obs_condition['1D'] and not obs_condition['VD'],
        '2D_scale': not obs_condition['1D'] and not any(obs_condition[key] for key in ['vector', 'SD', 'HD', 'VZ']),
        '3D_scale': obs_condition['3D'] and not any(obs_condition[key] for key in ['vector', 'SD'])
    }
    # Add constraints for translations
    for coord_label in controls.labels:
        if coord_label != 'o':
            translation_condition = np.zeros(controls.n_coords)
            translation_condition[coord_indices[coord_label]] = 1
            R.append(translation_condition)
            constraints_list.append(f'translation {coord_label}')
    # Add rotate condition
    if constraint_condition['rotate']:
        R.append(_rotate_condition(controls, coord_indices))
        constraints_list.append('rotate')
    # Add scale constraint
    if constraint_condition['1D_scale']:
        R.append(_scale_1D_condition(controls, coord_indices))
        constraints_list.append('1D scale')
    elif constraint_condition['2D_scale']:
        R.append(_scale_2D_condition(controls, coord_indices))
        constraints_list.append('2D scale')
    elif constraint_condition['3D_scale']:
        R.append(_scale_3D_condition(controls, coord_indices))
        constraints_list.append('3D scale')

    return {'R': np.vstack(R), 'constraints': constraints_list}


def calculate_sW(controls: Controls, calculate_sX=True,
                 default_pt_sigma: float = customizations.default_measurement_sigma['default']['sP']) -> dict:
    """
    Determines the sW and optionally sX matrices.

    Calculates the diagonal matrix of control points' weights (sW) and, optionally,
    the coefficient matrix for increments (sX) for adjustments that include control points' standard deviations.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - controls: (Controls): Coordinates of control points and their standard deviations.
    - calculate_sX: (bool): If True, the function calculates the sX matrix for coefficient matrix increments.
      Defaults to True.
    - default_pt_sigma: (float, optional): Default standard deviation of points location.
      Defaults as 'default sP' in customisations.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - dict: A dictionary containing matrices:
        - sW (np.ndarray): The weight matrix for the control points.
        - sX (np.ndarray, optional): The coefficient matrix for control points increments.
    """
    # Pre-compute reusable variables
    default_point_weight = len(controls.coord_labels) * inv_sq_error(default_pt_sigma)
    coord_indices = controls.coord_index
    # Initialize matrices
    sW = np.zeros((controls.n_coords, controls.n_coords))
    sX = np.eye(controls.n_coords) if calculate_sX else None
    # Iterate through control points and calculate weights
    for point, row in controls.iterrows():
        coord_idx = coord_indices[point]
        for coord, idx in coord_idx.items():
            s_coord = f's{coord}'
            coord_sigma = row.get(s_coord, np.nan)
            if s_coord == 'so':
                if sX is not None:
                    sX[idx, idx] = 0
                continue
            # Use provided sigma value if valid
            if not np.isnan(coord_sigma):
                if coord_sigma > 0:
                    sW[idx, idx] = inv_sq_error(coord_sigma)
                elif coord_sigma == -1:
                    sW[idx, idx] = 0
                continue
            sW[idx, idx] = default_point_weight
    # Save results to dictionary
    matrices = {'sW': sW}
    if sX is not None:
        matrices['sX'] = sX

    return matrices


def approx_orientation(controls: Controls, from_to_hz: pd.Series) -> pd.Series:
    """
    Computes the approximate orientation constant at control points.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - controls: (Controls): Coordinates of control points.
    - from_to_hz: (pd.Series): A Series with MultiIndex of ('FROM', 'TO') and values of horizontal directions.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - pd.Series: A Series containing the computed orientation constant for each point.
    """
    prev_hz_from = ''
    for (hz_from, hz_to), hz in from_to_hz.items():
        if prev_hz_from != hz_from:
            prev_hz_from = hz_from
            orientation = azimuth(controls.loc[hz_from], controls.loc[hz_to]) - hz
            if orientation < 0:
                orientation += 2 * np.pi
            controls.at[hz_from, 'o'] = orientation
    return controls['o']


def set_obs_weight(row: pd.Series, measurement_type: str, measurement_value: float,
                   default_sigma: dict = customizations.default_measurement_sigma['default']) -> float:
    """
    Computes the observation weight based on the measurement type and its standard deviation.

    If an estimated sigma value is given in the dataset, the point weight is calculated based on it.
    If this value is not present, the point weight is calculated based on the default value defined in
    the default_sigma dictionary.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - row: (pd.Series): A row from a measurements' dataset.
    - measurement_type: (str): The type of measurement.
    - measurement_value: (float): The value of the measurement.
    - default_sigma: (dict): Default standard deviations for measurements and points' localization.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - float: The observation weight.
    """
    sigma_label = f's{measurement_type}'
    measurement_sigma = row.get(sigma_label, np.nan)
    if not np.isnan(measurement_sigma) and measurement_sigma > 0:
        return inv_sq_error(measurement_sigma)
    elif sigma_label in Measurements.acceptable_labels['linear_measurements_sigma']:
        return inv_sq_error(default_sigma[sigma_label][0] + default_sigma[sigma_label][1] * measurement_value / 1000)
    return inv_sq_error(default_sigma[sigma_label])


def inv_sq_error(sigma: float) -> float:
    """
    Sets the observation weight based on the standard deviation.

    The weight is defined as the inverse of the square of the standard deviation.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - sigma: (float): Standard deviation value.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - float: Weight.
    """
    return np.divide(1, sigma**2)


# Helper functions *****************************************************************************************************
def _get_coords_differences(point_from: pd.Series, point_to: pd.Series) -> dict:
    """
    Calculates the coordinate differences between station and aim.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - point_from: (pd.Series): Coordinates of the station point.
    - point_to: (pd.Series): Coordinates of the aim point.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - dict: Coordinate differences between the points.
        - dx: (float, optional): The difference of the x-coordinates of the station and the aim point.
        - dy: (float, optional): The difference of the y-coordinates of the station and the aim point.
        - dz: (float, optional): The difference of the z-coordinates of the station and the aim point.
        - o: (float, optional): Value of orientation constant column at the station point.
    """
    differences = {f'd{label}': point_to[label] - point_from[label] for label in point_from.index}
    # Add orientation if is in points
    if 'o' in point_from.index:
        differences['o'] = point_from['o']

    return differences


def _rotate_condition(controls: Controls, coord_indices: dict) -> np.ndarray:
    """
    Creates a constraint for rotation around the center of the control network.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:
        controls (Controls): Coordinates of control points.
        coord_indices (dict): Dictionary with indexes of coordinate columns.

    --------------------------------------------------------------------------------------------------------------------
    Returns:
        np.ndarray: The rotation constraint.
    """
    rotate_cond = np.zeros(controls.n_coords)
    rotate_cond[coord_indices['x']] = controls['y'].values
    rotate_cond[coord_indices['y']] = -controls['x'].values
    return rotate_cond


def _scale_1D_condition(controls: Controls, coord_indices: dict) -> np.ndarray:
    """
    Creates a constraint for 1D scaling.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:
        controls (Controls): Coordinates of control points.
        coord_indices (dict): Dictionary with indexes of coordinate columns.

    --------------------------------------------------------------------------------------------------------------------
    Returns:
        np.ndarray: The 1D scale constraint.
    """
    scale = np.zeros(controls.n_coords)
    scale[coord_indices['z']] = controls['z'].values
    return scale


def _scale_2D_condition(controls: Controls, coord_indices: dict) -> np.ndarray:
    """
    Creates a constraint for 2D scaling.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:
        controls (Controls): Coordinates of control points.
        coord_indices (dict): Dictionary with indexes of coordinate columns.

    --------------------------------------------------------------------------------------------------------------------
    Returns:
        np.ndarray: The 2D scale constraint.
    """
    scale = np.zeros(controls.n_coords)
    scale[coord_indices['x']] = controls['x'].values
    scale[coord_indices['y']] = controls['y'].values
    return scale


def _scale_3D_condition(controls: Controls, coord_indices: dict) -> np.ndarray:
    """
    Creates a constraint for 3D scaling.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:
        controls (Controls): Coordinates of control points.
        coord_indices (dict): Dictionary with indexes of coordinate columns.

    --------------------------------------------------------------------------------------------------------------------
    Returns:
        np.ndarray: The 3D scale constraint.
    """
    scale = np.zeros(controls.n_coords)
    scale[coord_indices['x']] = controls['x'].values
    scale[coord_indices['y']] = controls['y'].values
    scale[coord_indices['z']] = controls['z'].values
    return scale
