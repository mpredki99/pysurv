# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import numpy as np
import pandas as pd
import pytest

from pysurv.modules.adjustment.matrices import equations_system
from pysurv.modules import Controls, Measurements

@pytest.fixture
def controls():
    return Controls(
            pd.DataFrame(
                {'ID': ['P01', 'P02', 'P03', 'P04'],
                 'x': [317.08, 706.43, 1063.29, 852.32],
                 'y': [1606.89, 1785.03, 1552.50, 1293.80],
                 'z': [103.57, 99.73, 99.19, 97.48]
                 })
    )

@pytest.fixture
def measurements():
    return Measurements(
            pd.DataFrame(
                {'FROM': ['P01', 'P01', 'P02', 'P02', 'P02', 'P03', 'P03', 'P04', 'P04', 'P04'],
                 'TO': ['P02', 'P03', 'P01', 'P03', 'P04', 'P01', 'P04', 'P01', 'P02', 'P03'],
                 'HD': [428.173, 748.195, np.nan, 425.933, 512.430, np.nan, 333.814, 620.088, np.nan, np.nan],
                 'sHD': [0.0024, 0.0027, np.nan, 0.0024, 0.0025, np.nan, 0.0023, 0.0026, np.nan, np.nan],
                 'HZ': [308.7758, 276.8270, 384.1063, 120.0253, 75.1682, 27.2014, 88.2790, 270.0645, 222.1380,
                        160.2063],
                 'VZ': [100.5715, 100.3728, 99.4291, 100.0820, 100.2799, 99.6268, 100.3253, 99.3742, 99.7206, 99.6756]
                 })
    )

def test_XY_keys(controls, measurements):
    methods = {
        'observations': 'ordinary',
        'obs_c': None,
        'free': 'ordinary',
        'free_c': None,
    }
    matrices = equations_system(controls, measurements, methods)
    assert list(matrices.keys()) == ['X', 'Y']
    assert all(isinstance(value, np.ndarray) for value in matrices.values())

def test_XYW_keys(controls, measurements):
    methods = {
        'observations': 'weighted',
        'obs_c': None,
        'free': 'ordinary',
        'free_c': None,
    }
    matrices = equations_system(controls, measurements, methods)
    assert list(matrices.keys()) == ['X', 'Y', 'W']
    assert all(isinstance(value, np.ndarray) for value in matrices.values())

def test_XY_sW_sX_keys(controls, measurements):
    methods = {
        'observations': 'ordinary',
        'obs_c': None,
        'free': None,
        'free_c': None,
    }
    matrices = equations_system(controls, measurements, methods)
    assert list(matrices.keys()) == ['X', 'Y', 'sW', 'sX']
    assert all(isinstance(value, np.ndarray) for value in matrices.values())

def test_XYW_sW_sX_keys(controls, measurements):
    methods = {
        'observations': 'weighted',
        'obs_c': None,
        'free': None,
        'free_c': None,
    }
    matrices = equations_system(controls, measurements, methods)
    assert list(matrices.keys()) == ['X', 'Y', 'W', 'sW', 'sX']
    assert all(isinstance(value, np.ndarray) for value in matrices.values())

def test_XY_R_sW_keys(controls, measurements):
    methods = {
        'observations': 'ordinary',
        'obs_c': None,
        'free': 'weighted',
        'free_c': None,
    }
    matrices = equations_system(controls, measurements, methods)
    assert list(matrices.keys()) == ['X', 'Y', 'R', 'constraints', 'sW']
    assert all(isinstance(matrices[value], np.ndarray) for value in matrices if value != 'constraints')
    assert isinstance(matrices['constraints'], list)

def test_XYW_R_sW_keys(controls, measurements):
    methods = {
        'observations': 'weighted',
        'obs_c': None,
        'free': 'weighted',
        'free_c': None,
    }
    matrices = equations_system(controls, measurements, methods)
    assert list(matrices.keys()) == ['X', 'Y', 'W', 'R', 'constraints', 'sW']
    assert all(isinstance(matrices[value], np.ndarray) for value in matrices if value != 'constraints')
    assert isinstance(matrices['constraints'], list)

def test_robust_keys(controls, measurements):
    methods = {
        'observations': 'huber',
        'obs_c': None,
        'free': 'tukey',
        'free_c': None,
    }
    matrices = equations_system(controls, measurements, methods)
    assert list(matrices.keys()) == ['X', 'Y', 'W', 'R', 'constraints', 'sW']
    assert all(isinstance(matrices[value], np.ndarray) for value in matrices if value != 'constraints')
    assert isinstance(matrices['constraints'], list)
