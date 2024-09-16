# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
This module defines default configurations for various parameters used in the adjustment of surveying control networks.

The following dictionaries are provided:

1. `methods`: Specifies default methods for adjusting observations and reference points in free adjustment.
    - `observations`: Default method for adjusting observations (str).
    - `obs_c`: Tuning constant for robust adjustment of observations (float, tuple of two floats, or None).
    - `free`: Default method for free adjustment (str or None).
    - `free_c`: Tuning constant for robust free adjustment (float, tuple of two floats, or None).

    If a robust method is provided and the tuning constant is None,
    a default theoretical value of the method will be used to calculate reweighting factors.

2. `default_measurement_sigma`: Contains default measurement sigma values for various quantities.
    - `sSD`: Default sigma value of 3D spatial distances (tuple of two floats).
    - `sHD`: Default sigma value of 2D horizontal distances (tuple of two floats).
    - `sVD`: Default sigma value of 1D vertical distances (float).
    - `sdx`: Default sigma value of x-component of GNSS vectors (tuple of two floats).
    - `sdy`: Default sigma value of y-component of GNSS vectors (tuple of two floats).
    - `sdz`: Default sigma value of z-component of GNSS vectors (tuple of two floats).
    - `sA`: Default sigma value of azimuthal angles (float).
    - `sHZ`: Default sigma value of horizontal directions (float).
    - `sVZ`: Default sigma value of vertical zenith angles (float).
    - `sVH`: Default sigma value of vertical horizontal angles (float).
    - `sP`: Default sigma value of the position of control points (float).

3. `iterate_params`: Defines default parameters for iteration process in adjustment computations.
    - `tolerance`: Convergence tolerance (float).
    - `max_iter`: Maximum number of iterations (int).

4. `report_params`: Defines parameters used to create report with adjustment results.
    - `approx_precision`: number of decimals of approx coordinates and linear measurements (int).
    - `adjusted_precision`: number of decimals of adjusted coordinates and linear measurements (int).
    - `angle_unit`: unit to represent angles (str).
    - `angle_precision`: number of decimals of angles (int).
"""
methods = dict(
    default={'observations': 'weighted',
             'obs_c': None,
             'free': None,
             'free_c': None,
             },
)
default_measurement_sigma = dict(
    default={'sSD': (0.003, 0.002),
             'sHD': (0.003, 0.002),
             'sVD': (0.003, 0.002),
             'sdx': (0.002, 0.001),
             'sdy': (0.002, 0.001),
             'sdz': (0.002, 0.001),
             'sHZ': 0.0020,
             'sVZ': 0.0020,
             'sVH': 0.0020,
             'sA': 0.0020,
             'sP': 0.05,
             },
)
iterate_params = dict(
    default={'tolerance': 0.00001,
             'max_iter': 100,
             },
)
report_params = dict(
    default={'approx_precision': 3,
             'adjusted_precision': 4,
             'angle_unit': 'grad',
             'angle_precision': 4,
             },
)
