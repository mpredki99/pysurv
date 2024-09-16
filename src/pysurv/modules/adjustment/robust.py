# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
This module provides various weight functions used in robust estimation.

The module includes four main categories of weight functions:
1. With Tolerance: Functions that do not change the weights of observations that fall within a specified range.
2. Bell Curves: Functions based on different types of bell-shaped curves.
3. Trigonometric Functions: Functions that incorporate trigonometric functions, including hyperbolic functions.
4. Others: Miscellaneous functions that do not fall into the previous categories but provide unique methods for adjusting weights.

Each function takes normalized residuals as input and returns an array of reweight factors.
Theoretical values of tuning constants are provided as default.
"""
from typing import Tuple
import numpy as np
from scipy.special import erf

# Theoretical tuning constants for various weight functions
theoretical_c = {"huber": 1.345,
                 "slope": (2, 2),
                 "hampel": (1.7, 3.4, 8.5),
                 "danish": 2.5,
                 "epanechnikov": (3.674, 2),
                 "tukey": (4.685, 2),
                 "jacobi": (4.687, 1),
                 "exponential": (2, 2),
                 "sigma": (2, 2),
                 "error_func": (1.414, 2),
                 "cauchy": (2.385, 2),
                 "t": (1, 2),
                 "chain_bell": (1, 1),
                 "chain": 1,
                 "andrews": 4.207,
                 "wave": 2.5,
                 "half_wave": 2.5,
                 "wigner": 3.137,
                 "ellipse_curve": 2.5,
                 "trim": 2.5,
                 }


# WITH TOLERANCE
def huber(v: np.ndarray, c: float = theoretical_c['huber']) -> np.ndarray:
    """
    Huber M-estimator weight function.

    Re-weighting coefficients for residuals are calculated as:
    1       : for: v <= c
    c / |v| : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 1.345.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    return np.where(v > c, c / v, 1)


def slope(v: np.ndarray, c: Tuple[float, float] = theoretical_c['slope']) -> np.ndarray:
    """
    Slope weight function.

    Re-weighting coefficients for residuals are calculated as:
    1                 : for: v <= c
    1 + (c - |v|) / a : for: v > c and coefficients greater than 0,
    the rest of them take the value 0.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of two floats): Tuning constants (c, a), default is (2, 2).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, a = c
    weights = 1 + (c - v) / a
    return np.clip(weights, 0, 1)


def hampel(v: np.array, c: Tuple[float, float, float] = theoretical_c['hampel']) -> np.ndarray:
    """
    Hampel weight function.

    Re-weighting coefficients for residuals are calculated as:
    1                             : for: v <= a
    a / |v|                       : for: a < v <= b
    a / |v| * (c - |v|) / (c - b) : for: b < v <= c
    0                             : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of three floats): Tuning constants (a, b, c), default is (1.7, 3.4, 8.5).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    a, b, c = c
    weights = np.ones_like(v)
    weights[v > a] = np.divide(a, v[v > a])
    weights[v > b] = np.divide(a, v[v > b]) * np.divide(c - v[v > b], c - b)
    weights[v > c] = 0
    return weights


def danish(v: np.ndarray, c: float = theoretical_c['danish']) -> np.ndarray:
    """
    Danish weight function.

    Re-weighting coefficients for residuals are calculated as:
    1           : for: v <= c
    exp(-v / c) : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 2.5.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    return np.where(v > c, np.exp(-v / c), 1)


# BELL CURVES
def epanechnikov(v: np.ndarray, c: Tuple[float, float] = theoretical_c['epanechnikov']) -> np.ndarray:
    """
    Epanechnikov weight function.

    Re-weighting coefficients for residuals are calculated as:
    1 - (v / c)^n : for: v <= c
    0             : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of two floats): Tuning constants (c, n), default is (3.674, 2.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

      - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return np.where(v <= c, 1 - (v / c)**n, 0)


def tukey(v: np.ndarray, c: Tuple[float, float] = theoretical_c['tukey']) -> np.ndarray:
    """
    Tukey weight function.

    Re-weighting coefficients for residuals are calculated as:
    (1 - (v / c)^n)^n : for: v <= c
    0                 : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of two floats): Tuning constants (c, n), default is (4.685, 2.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return np.where(v <= c, (1 - (v / c)**n)**n, 0)


def jacobi(v: np.ndarray, c: Tuple[float, float] = theoretical_c['jacobi']) -> np.ndarray:
    """
    Jacobi weight function.

    Re-weighting coefficients for residuals are calculated as:
    (1 - (v / c)^n)^n * (1 + (v / c)^n)^n : for: v <= c
    0                                     : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of two floats): Tuning constants (c, n), default is (4.687, 1.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return np.where(v <= c, (1 - (v / c)**n)**n * (1 + (v / c)**n)**n, 0)


def exponential(v: np.ndarray, c: Tuple[float, float] = theoretical_c['exponential']) -> np.ndarray:
    """
    Exponential weight function.

    Re-weighting coefficients for residuals are calculated as:
    exp((-v / c)^n)

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of two floats): Tuning constants (c, n), default is (2.0, 2.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return np.exp(-(v / c)**n)


def sigma(v: np.ndarray, sigma_sq: float, c: Tuple[float, float] = theoretical_c['sigma']) -> np.ndarray:
    """
    Sigma weight function.

    Re-weighting coefficients for residuals are calculated as:
    exp(-v^n / (sigma_sq * c))

    The shape of the function changes with each iteration as the value of the residual variance sigma_sq changes.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - sigma_sq: (float): Residual variance value.
    - c: (tuple of two floats): Tuning constants (c, n), default is (2.0, 2.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return np.exp(-v**n / (sigma_sq * c))


def error_func(v: np.ndarray, c: Tuple[float, float] = theoretical_c['error_func']) -> np.ndarray:
    """
    Error function weight function.

    Re-weighting coefficients for residuals are calculated as:
    1 - erf((v / c)^n)

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of two floats): Tuning constants (c, n), default is (1.414, 2.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return 1 - erf((v / c)**n)


def cauchy(v: np.ndarray, c: Tuple[float, float] = theoretical_c['cauchy']) -> np.ndarray:
    """
    Cauchy weight function.

    Re-weighting coefficients for residuals are calculated as:
    1 / (1 + (v / c)^n)

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of two floats): Tuning constants (c, n), default is (2.385, 2.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return np.divide(1, 1 + (v / c)**n)


def t(v: np.ndarray, k: int, c: Tuple[float, float] = theoretical_c['t']) -> np.ndarray:
    """
    Student's t weight function.

    Re-weighting coefficients for residuals are calculated as:
    (1 + (v^n) / (c  k))^(-(k + 1) / 2)

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - k: (int): Degrees of freedom.
    - c: (tuple of two floats): Tuning constants (c, n), default is (1.0, 2.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return np.power(1 + v**n / (c * k), -(k + 1) / 2)


def chain_bell(v: np.ndarray, c: Tuple[float, float] = theoretical_c['chain_bell']) -> np.ndarray:
    """
    Chain Bell weight function.

    Re-weighting coefficients for residuals are calculated as:
    1 / (cosh((v^n * e) / (2 * c)))

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (tuple of two floats): Tuning constants (c, n), default is (1.0, 1.0).

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c, n = c
    return np.divide(1, np.cosh(np.divide(v**n * np.e, 2 * c)))


# TRIGONOMETRIC FUNCTIONS
def chain(v: np.ndarray, c: float = theoretical_c['chain']) -> np.ndarray:
    """
    Chain weight function.

    Re-weighting coefficients for residuals are calculated as:
    -cosh((v * e) / (2 * c)) + 2 : for coefficients greater than 0,

    the rest of them take the value 0.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 1.0.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    weights = -np.cosh((v * np.e) / (2 * c)) + 2
    return np.where(weights < 0, 0, weights)


def andrews(v: np.ndarray, c: float = theoretical_c['andrews']) -> np.ndarray:
    """
    Andrews weight function.

    Re-weighting coefficients for residuals are calculated as:
    sinc(v / c) : for: v <= c
    0           : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 4.207.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    return np.where(v <= c, np.sinc(v / c), 0)


def wave(v: np.ndarray, c: float = theoretical_c['wave']) -> np.ndarray:
    """
    Wave weight function.

    Re-weighting coefficients for residuals are calculated as:
    cos((v * π) / c) + 1) / (2) : for: v <= c
    0                           : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 2.5.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    return np.where(v <= c, np.divide(np.cos(v * np.pi / c) + 1, 2), 0)


def half_wave(v: np.ndarray, c: float = theoretical_c['half_wave']) -> np.ndarray:
    """
    Half Wave weight function.

    Re-weighting coefficients for residuals are calculated as:
    cos((v * pi) / (2  c)) : for: v <= c
    0                      : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 2.5.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    return np.where(v <= c, np.cos(v * np.pi / (2 * c)), 0)


# OTHERS
def wigner(v: np.ndarray, c: float = theoretical_c['wigner']) -> np.ndarray:
    """
    Wigner weight function.

    Re-weighting coefficients for residuals are calculated as:
    sqrt(1 - (v / c)^2) : for: v <= c
    0                   : for: v > c

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 3.137.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    return np.where(v <= c, np.sqrt(1 - (v / c)**2), 0)


def ellipse_curve(v: np.ndarray, c: float = 2.5) -> np.ndarray:
    """
    Ellipse Curve weight function.

    Re-weighting coefficients for residuals are calculated as:
    (1 + sqrt(1 - (v / c)^2)) / 2       : for: v <= c
    (1 - sqrt(1 - ((v - 2c) / c)^2) / 2 : for: c > v <= 2c
    0                                   : for: v > 2c

    ---------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 2.5.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    c2 = 2 * c
    weights = np.ones_like(v)
    mask1 = v <= c
    mask2 = (v > c) & (v <= c2)
    weights[mask1] = np.divide(1 + np.sqrt(1 - (v[mask1] / c) ** 2), 2)
    weights[mask2] = np.divide(1 - np.sqrt(1 - ((v[mask2] - c2) / c) ** 2), 2)
    weights[v > c2] = 0
    return weights


def trim(v: np.ndarray, c: float = theoretical_c['trim']) -> np.ndarray:
    """
    Trim weight function.

    Rejects from the set of observations, those whose values of normalized residuals
    exceed the value of the tuning constant c.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - v: (np.ndarray): Normalized values of residuals.
    - c: (float): Tuning constant, default is 2.5.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - np.ndarray: Parameters for adjusting the weights.
    """
    v = np.abs(v)
    return np.where(v <= c, 1, 0)


# Dict of module's weight functions
methods = {"huber": huber,
           "slope": slope,
           "hampel": hampel,
           "danish": danish,
           "epanechnikov": epanechnikov,
           "tukey": tukey,
           "jacobi": jacobi,
           "exponential": exponential,
           "sigma": sigma,
           "error_func": error_func,
           "cauchy": cauchy,
           "t": t,
           "chain_bell": chain_bell,
           "chain": chain,
           "andrews": andrews,
           "wave": wave,
           "half_wave": half_wave,
           "wigner": wigner,
           "ellipse_curve": ellipse_curve,
           "trim": trim,
           }
