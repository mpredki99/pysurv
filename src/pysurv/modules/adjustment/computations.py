# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
This module contains functions to perform least squares adjustments for surveying control networks.

The adjustment process involves iterating over observation equations, updating weights
based on robust methods, and solving task with free adjustment approach.
The results include the adjusted coordinates, covariance matrices and other information
about adjustment process.
"""
import warnings
from typing import Tuple
import numpy as np
from numpy.linalg import pinv
from pysurv.modules import Controls, Measurements
from pysurv import customizations
from . import robust
from .matrices import equations_system


def adjust(controls: Controls, measurements: Measurements, matrices: dict,
           methods: dict = customizations.methods['default'],
           iterate_params: dict = customizations.iterate_params['default']) -> dict:
    """
    Perform least squares adjustment on a set of control points and measurements.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - controls: (Controls): Coordinates of control points.
    - measurements: (Measurements): The measurement dataset.
    - matrices: (dict): A dictionary containing matrices X, Y, and optionally W, R, sW, sX.
    - methods: (dict): Dictionary with methods to use to set observation and control point weights,
      and also with tuning constants. Defaults as 'default' in customisations.
    - iterate_params: (dict): Parameters for iteration process in computations. Defaults as 'default' in customisations.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - dict: A dictionary containing the adjustment results:
        - n_iter: (int): The number of iterations performed.
        - n_measurements: (int): Number of measurements.
        - n_fixed_coords: (int): Number of fixed reference coordinates.
        - n_sigma_coords: (int): Number of movable reference coordinates.
        - n_unknowns: (int): Number of unknowns.
        - r_norm: (int): Normalized residuals.
        - b_norm: (int): Normalized increments.
        - sigma_zero: (list): List of sigma zero values from each iteration.
        - pt_sigma_zero: (list): List of control points increments sigma zero values from each iteration.
        - approx_coordinates: (Controls): The approximate control points coordinates before adjustment.
        - adjusted_coordinates: (Controls): The adjusted control points coordinates.
        - constraints_list: (list): List of inner constraints applied to the network.
    """
    # Get values of matrices
    X, R, sW, sX = matrices.get('X'), matrices.get('R'), matrices.get('sW'), matrices.get('sX')
    # Determine degrees of freedom
    n_constraints = len(R) if R is not None else 0
    n_measurements = len(X)
    n_fixed_coords = ((controls[controls.sigma_labels] < iterate_params['tolerance']) &
                      (controls[controls.sigma_labels] != -1)).sum().sum().astype(int)
    n_sigma_coords = len(sW.diagonal()[sW.diagonal() != 0]) - n_fixed_coords if sW is not None else (
                     controls[controls.coord_labels].notna().count().sum().astype(int))
    n_unknowns = X.shape[1] - n_fixed_coords
    # Add information about degrees of freedom to matrices dictionary
    matrices['k'] = int(n_measurements + n_constraints - n_unknowns)
    # Store coordinates from controls file
    approx_coordinates = controls.copy_with_type()
    # Initialize reusable values
    sigma_zero_sq, pt_sigma_zero_sq = list(), list()
    sum_increments = np.zeros([X.shape[1], 1])
    threshold = -np.log10(min(iterate_params['tolerance'],
                              customizations.iterate_params['default']['tolerance'])**2).astype(int) + 1
    n_iter = 0
    while True:
        n_iter += 1
        # Update the weights based on robust methods
        if n_iter > 1:
            update_weights(matrices, adjustment_results, methods, pt_sigma_zero_sq[-1])
        # Solve system of observation equations
        try:
            adjustment_results = iterate(matrices)
        except np.linalg.LinAlgError:
            warnings.warn(f'SVD did not converge in {n_iter} iteration. The calculations have been aborted!')
            _finalize_results(adjustment_results, n_iter, n_measurements, n_fixed_coords, n_sigma_coords, n_unknowns,
                sigma_zero_sq, pt_sigma_zero_sq, approx_coordinates, controls, matrices.get('constraints'))
            return adjustment_results
        # Continue calculations if SVD converged
        increments = adjustment_results['increments']
        sum_increments += increments
        # Calculate references sigma zero
        sigma_zero_sq.append(adjustment_results['sigma_zero_sq'])
        if sW is not None and n_sigma_coords > 0:
            pt_sigma_zero_sq.append(np.divide((sum_increments.reshape(-1)**2 * sW.diagonal()).sum(), n_sigma_coords))
        elif sW is None and n_sigma_coords > 0:
            pt_sigma_zero_sq.append(np.divide((sum_increments**2).sum(), n_sigma_coords))
        else:
            pt_sigma_zero_sq.append(1)
        # Apply threshold
        adjustment_results['cov_b'] = np.round(adjustment_results['cov_b'], decimals=threshold)
        adjustment_results['cov_Y'] = np.round(adjustment_results['cov_Y'], decimals=threshold)
        adjustment_results['cov_r'] = np.round(adjustment_results['cov_r'], decimals=threshold)
        # Update coordinates unless it is first iteration of robust free adjustment
        if methods['free'] not in robust.methods or n_iter != 1:
            # Get the row and column indices of the non-NaN values
            row_indices, col_indices = np.where(~controls[controls.labels].isna())
            # Apply the increments to the corresponding positions
            for row, col, increment in zip(row_indices, col_indices, increments):
                controls.loc[controls.index[row], controls.labels[col]] += increment
        # Check for convergence based on iteration parameters
        if (
                n_iter >= iterate_params['max_iter'] or
                np.all(adjustment_results['increments'] <= iterate_params['tolerance'])
        ):
            _finalize_results(
                adjustment_results,
                n_iter,
                n_measurements,
                n_fixed_coords,
                n_sigma_coords,
                n_unknowns,
                sigma_zero_sq,
                pt_sigma_zero_sq,
                approx_coordinates,
                controls,
                matrices.get('constraints')
            )
            return adjustment_results
        # Update matrices for the next iteration
        new_XY = equations_system(
            controls,
            measurements,
            methods={
                'observations': 'ordinary',
                'obs_c': None,
                'free': 'ordinary',
                'free_c': None
            }
        )
        matrices.update({'X': new_XY['X'], 'Y': new_XY['Y']})


def iterate(matrices: dict) -> dict:
    """
    Perform iteration for least squares adjustment.

    Solves a system of equations and determines covariance matrices.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - matrices: (dict): A dictionary containing matrices:
        - X: (np.ndarray): The coefficient matrix.
        - Y: (np.ndarray): The vector of observed values.
        - W: (np.ndarray, optional): The observations weight matrix.
        - R: (np.ndarray, optional): The matrix of network inner constraints.
        - sX: (np.ndarray, optional): The coefficient matrix for control points increments.
        - sW: (np.ndarray, optional): The weight matrix for the control points.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - dict: A dictionary containing the iteration results:
        - increments: (np.ndarray): Vector of calculated increments to coordinates.
        - obs_residuals: (np.ndarray): Vector of observation residuals.
        - sigma_zero_sq: (np.ndarray): Residual variance.
        - N: (np.ndarray): Matrix of normal equations.
        - L: (np.ndarray): Vector of dependent variable.
        - cov_b: (np.ndarray): Covariance matrix of increments.
        - cov_Y: (np.ndarray): Covariance matrix of observations.
        - cov_r: (np.ndarray): Covariance matrix of residuals.
    """
    # Get values of matrices
    X, W, Y = matrices.get('X'), matrices.get('W'), matrices.get('Y')
    sW, R, sX = matrices.get('sW'), matrices.get('R'), matrices.get('sX')
    k = matrices.get('k')
    # Calculate normal equations
    N, L = normal_equations(X, W, Y, sW, R, sX)
    # Cofactor matrix
    Qbb = pinv(N)
    # Solve system
    b = Qbb @ L
    # Calculate residuals for observations
    r = X @ b - Y
    # Calculate residual variance
    if W is not None and k > 0:
        sigma_zero_sq = np.divide((r.reshape(-1)**2 * W.diagonal()).sum(), k)
    elif W is None and k > 0:
        sigma_zero_sq = np.divide((r**2).sum(), k)
    else:
        sigma_zero_sq = 1
    # Determine covariance matrices
    cov_b = sigma_zero_sq * Qbb
    cov_Y = X @ cov_b @ X.T
    cov_r = (sigma_zero_sq * (pinv(W) if W is not None else np.eye(len(cov_Y)))) - cov_Y

    return {'increments': b,
            'obs_residuals': r,
            'sigma_zero_sq': sigma_zero_sq,
            'N': N,
            'L': L,
            'cov_b': cov_b,
            'cov_Y': cov_Y,
            'cov_r': cov_r,
            }


def normal_equations(X: np.ndarray, W: np.ndarray, Y: np.ndarray,
                     sW: np.ndarray, R: np.ndarray, sX: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Determine the normal equations and dependent variable for the adjustment process.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - X: (np.ndarray): Coefficient matrix.
    - W: (np.ndarray): The observations weight matrix.
    - Y: (np.ndarray): The vector of observed values.
    - sW: (np.ndarray): The weight matrix for the control points.
    - R: (np.ndarray): The matrix of network inner constraints.
    - sX: (np.ndarray): The coefficient matrix for control points increments.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - Tuple containing normal equations and dependent variable:
        - np.ndarray: The matrix of normal equations.
        - np.ndarray: The vector of dependent variable.
    """
    # Calculate first factor of normal equations
    N1 = X.T @ (W @ X if W is not None else X)
    # Calculate second factor of normal equations
    if sX is not None and sW is not None:
        N2 = sX.T @ sW @ sX
    elif R is not None and sW is not None:
        N2 = R.T @ R @ sW
    else:
        N2 = None
    # Calculate normal equations
    N = N1 + N2 if N2 is not None else N1
    L = X.T @ (W @ Y if W is not None else Y)
    return N, L


def update_weights(matrices: dict, iterate_results: dict,
                   methods: dict, pt_sigma_zero_sq: float) -> None:
    """
    Update the weight matrices based on robust methods during the iteration process.

    This function modifies the input array in-place without returning a value.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - matrices: (dict): A dictionary containing weight matrices.
    - iterate_results: (dict): The results from the current iteration.
    - methods: (dict): Dictionary specifying the robust methods and tuning constants to use.
    - pt_sigma_zero_sq: (float): Control points' increments squared sigma zero parameter.
    """
    if methods['observations'] in robust.methods:
        kwargs = dict()
        # Get user's provided tuning constant or don't create to get default value
        if methods['obs_c'] is not None:
            kwargs['c'] = methods['obs_c']
        # Add another keyword arguments if necessary
        if methods['observations'] == 'sigma':
            kwargs['sigma_sq'] = iterate_results.get('sigma_zero_sq')
        elif methods['observations'] == 't':
            kwargs['k'] = matrices.get('k')

        reweight(
            matrices.get('W'),
            iterate_results.get('obs_residuals').reshape(-1),
            iterate_results.get('cov_r').diagonal(),
            methods['observations'],
            kwargs
        )

    if methods['free'] in robust.methods:
        kwargs = dict()
        # Get user's provided tuning constant or don't create to get default value
        if methods['free_c'] is not None:
            kwargs['c'] = methods['free_c']
        # Add another keyword arguments if necessary
        if methods['free'] == 'sigma':
            kwargs['sigma_sq'] = pt_sigma_zero_sq
        elif methods['free'] == 't':
            kwargs['k'] = matrices.get('k')

        reweight(
            matrices.get('sW'),
            iterate_results.get('increments').reshape(-1),
            iterate_results.get('cov_b').diagonal(),
            methods['free'],
            kwargs
        )


def reweight(weights: np.ndarray, residuals: np.ndarray,
             residuals_variances: np.ndarray, method: str, kwargs: dict) -> None:
    """
    Recalculate the weights based on the normalized residuals and a selected robust method.

    This function modifies the input array in-place without returning a value.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - weights: (np.ndarray): The current weight matrix.
    - residuals: (np.ndarray): The residual values to be reweighted.
    - residuals_variances: (np.ndarray): The variances of the residuals.
    - method: (str): The robust method to use for reweighting.
    - **kwargs: (dict): Keyword arguments required by the robust method.
    """
    norm_residuals = normalize_residuals(residuals, residuals_variances)
    np.fill_diagonal(weights, weights.diagonal() * robust.methods[method](norm_residuals, **kwargs))


def normalize_residuals(residuals: np.ndarray, residuals_variances: np.ndarray) -> np.ndarray:
    """
    Normalize values of residuals based on their values and variances.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - residuals: (np.ndarray): The residual values.
    - residuals_variances: (np.ndarray): The variances of the residuals.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Array of normalized values of residuals.
    """
    invalid_variances = residuals_variances[residuals_variances < 0]
    if invalid_variances.size > 0:
        warnings.warn(f'Invalid variances encountered in sqrt: {invalid_variances}')
    with np.errstate(divide='ignore', invalid='ignore'):
        residuals_sigma = np.sqrt(residuals_variances)
        return np.where(residuals_sigma == 0 , np.nan, residuals / residuals_sigma)


# Helper functions *****************************************************************************************************
def _finalize_results(adjustment_results: dict,
                      n_iter: int,
                      n_measurements: int,
                      n_fixed_coords: int,
                      n_sigma_coords: int,
                      n_unknowns: int,
                      sigma_zero: list,
                      pt_sigma_zero: list,
                      approx_coordinates: Controls,
                      adjusted_coordinates: Controls,
                      constraints_list: list) -> None:
    """
    Add information to the results after the adjustment process, normalizing residuals and increments.

    Information are added in original object in-place without returning any value.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - adjustment_results: (dict): The results from the adjustment computations.
    - n_iter: (int): The number of iterations performed.
    - n_measurements: (int): Number of measurements.
    - n_fixed_coords: (int): Number of fixed reference coordinates.
    - n_sigma_coords: (int): Number of movable reference coordinates.
    - n_unknowns: (int): Number of unknowns.
    - r_norm: (int): Normalized residuals.
    - b_norm: (int): Normalized increments.
    - sigma_zero: (list): List of sigma zero values from each iteration.
    - pt_sigma_zero: (list): List of control points increments sigma zero values from each iteration.
    - approx_coordinates: (Controls): The approximate control points coordinates before adjustment.
    - adjusted_coordinates: (Controls): The adjusted control points coordinates.
    - constraints_list: (list, optional): List of inner constraints applied to the network.
    """
    # Normalize residuals
    r = adjustment_results.get('obs_residuals')
    cov_r = adjustment_results.get('cov_r')
    r_norm = normalize_residuals(r.reshape(-1), np.diagonal(cov_r))
    # Normalize increments
    b = adjustment_results.get('increments')
    cov_b = adjustment_results.get('cov_b')
    b_norm = normalize_residuals(b.reshape(-1), np.diagonal(cov_b))
    # Add info to the results
    adjustment_results['n_iter'] = n_iter
    adjustment_results['n_measurements'] = n_measurements
    adjustment_results['n_fixed_coords'] = n_fixed_coords
    adjustment_results['n_sigma_coords'] = n_sigma_coords
    adjustment_results['n_unknowns'] = n_unknowns
    adjustment_results['r_norm'] = r_norm
    adjustment_results['b_norm'] = b_norm
    adjustment_results['sigma_zero'] = np.sqrt(np.array(sigma_zero))
    adjustment_results['pt_sigma_zero'] = np.sqrt(np.array(pt_sigma_zero))
    adjustment_results['approx'] = approx_coordinates
    adjustment_results['adjusted'] = adjusted_coordinates
    if constraints_list is not None:
        adjustment_results['constraints'] = constraints_list
