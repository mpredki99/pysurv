# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import numpy as np
import pandas as pd
import pytest
from pysurv.modules import Measurements
from pysurv.modules.basic import to_rad

@pytest.fixture
def standard_data():
    return pd.DataFrame(
            {'FROM': ['P01', 'P01', 'P02', 'P02', 'P02', 'P03', 'P03', 'P04', 'P04', 'P04'],
             'TO': ['P02', 'P03', 'P01', 'P03', 'P04', 'P01', 'P04', 'P01', 'P02', 'P03'],
             'HD': [428.173, 748.195, np.nan, 425.933, 512.430, np.nan, 333.814, 620.088, np.nan, np.nan],
             'sHD': [0.0024, 0.0027, np.nan, 0.0024, 0.0025, np.nan, 0.0023, 0.0026, np.nan, np.nan],
             'HZ': [277.8982, 249.1443, 345.6957, 108.0228, 67.6514, 24.4813, 79.4511, 243.0581, 199.9242, 144.1857],
             'sHZ': [0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018],
             'VZ': [90.5143, 90.3355, 89.4862, 90.0738, 90.2519, 89.6641, 90.2928, 89.4368, 89.7485, 89.7080]
             })

@pytest.fixture
def non_standard_data():
    return pd.DataFrame(
            {'from': ['P01', 'P01', 'P02', 'P02', 'P02', 'P03', 'P03', 'P04', 'P04', 'P04'],
             'to': ['P02', 'P03', 'P01', 'P03', 'P04', 'P01', 'P04', 'P01', 'P02', 'P03'],
             'hd': [428.173, 748.195, np.nan, 425.933, 512.430, np.nan, 333.814, 620.088, np.nan, np.nan],
             'shd': [0.0024, 0.0027, np.nan, 0.0024, 0.0025, np.nan, 0.0023, 0.0026, np.nan, np.nan],
             'hz': [277.8982, 249.1443, 345.6957, 108.0228, 67.6514, 24.4813, 79.4511, 243.0581, 199.9242, 144.1857],
             'shz': [0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018],
             'vz': [90.5143, 90.3355, 89.4862, 90.0738, 90.2519, 89.6641, 90.2928, 89.4368, 89.7485, 89.7080]
             })

def test_initialization_standard_columns(standard_data):
    measurements = Measurements(standard_data)
    pd.testing.assert_index_equal(measurements.columns, pd.Index(['HD', 'sHD', 'HZ', 'sHZ', 'VZ']))

def test_initialization_with_nonstandard_columns(non_standard_data):
    measurements = Measurements(non_standard_data)
    pd.testing.assert_index_equal(measurements.columns, pd.Index(['HD', 'sHD', 'HZ', 'sHZ', 'VZ']))

def test_angle_conversion_to_radians(non_standard_data):
    measurements = Measurements(non_standard_data, angle_unit='deg')
    assert pytest.approx(measurements['HZ'][0]) == to_rad(277.8982, unit='deg')
    assert pytest.approx(measurements['VZ'][0]) == to_rad(90.5143, unit='deg')

def test_invalid_angle_unit(standard_data):
    with pytest.raises(ValueError):
        Measurements(standard_data, angle_unit='invalid_unit')

def test_types_property(non_standard_data):
    measurements = Measurements(non_standard_data)
    assert list(measurements.types) == ['HD', 'HZ', 'VZ']

def test_sigma_property(non_standard_data):
    measurements = Measurements(non_standard_data)
    assert list(measurements.sigma) == ['sHD', 'sHZ']

def test_linear_property(non_standard_data):
    measurements = Measurements(non_standard_data)
    assert list(measurements.linear) == ['HD']

def test_linear_sigma_property(non_standard_data):
    measurements = Measurements(non_standard_data)
    assert list(measurements.linear_sigma) == ['sHD']

def test_angular_property(non_standard_data):
    measurements = Measurements(non_standard_data)
    assert list(measurements.angular) == ['HZ', 'VZ']

def test_angular_sigma_property(non_standard_data):
    measurements = Measurements(non_standard_data)
    assert list(measurements.angular_sigma) == ['sHZ']

def test_invalid_column_names():
    invalid_data = pd.DataFrame({
        'from': ['P1', 'P2'],
        'to': ['P2', 'P1'],
        'dist': [10.0, 15.0],
        'ang': [20.0, 25.0],
        'time': [5.0, 7.0],
    })
    measurements = Measurements(invalid_data)
    assert measurements.empty

def test_itermeasurements(non_standard_data):
    measurements = Measurements(non_standard_data, angle_unit='deg')
    measurement_list = list(measurements.itermeasurements())
    expected_measurements = [
        (('P01', 'P02'), 'HD', 428.173),
        (('P01', 'P02'), 'HZ', to_rad(277.8982, unit='deg')),
        (('P01', 'P02'), 'VZ', to_rad(90.5143, unit='deg')),
        (('P01', 'P03'), 'HD', 748.195),
        (('P01', 'P03'), 'HZ', to_rad(249.1443, unit='deg')),
        (('P01', 'P03'), 'VZ', to_rad(90.3355, unit='deg')),
        (('P02', 'P01'), 'HZ', to_rad(345.6957, unit='deg')),
        (('P02', 'P01'), 'VZ', to_rad(89.4862, unit='deg')),
        (('P02', 'P03'), 'HD', 425.933),
        (('P02', 'P03'), 'HZ', to_rad(108.0228, unit='deg')),
        (('P02', 'P03'), 'VZ', to_rad(90.0738, unit='deg')),
        (('P02', 'P04'), 'HD', 512.43),
        (('P02', 'P04'), 'HZ', to_rad(67.6514, unit='deg')),
        (('P02', 'P04'), 'VZ', to_rad(90.2519, unit='deg')),
        (('P03', 'P01'), 'HZ', to_rad(24.4813, unit='deg')),
        (('P03', 'P01'), 'VZ', to_rad(89.6641, unit='deg')),
        (('P03', 'P04'), 'HD', 333.814),
        (('P03', 'P04'), 'HZ', to_rad(79.4511, unit='deg')),
        (('P03', 'P04'), 'VZ', to_rad(90.2928, unit='deg')),
        (('P04', 'P01'), 'HD', 620.088),
        (('P04', 'P01'), 'HZ', to_rad(243.0581, unit='deg')),
        (('P04', 'P01'), 'VZ', to_rad(89.4368, unit='deg')),
        (('P04', 'P02'), 'HZ', to_rad(199.9242, unit='deg')),
        (('P04', 'P02'), 'VZ', to_rad(89.7485, unit='deg')),
        (('P04', 'P03'), 'HZ', to_rad(144.1857, unit='deg')),
        (('P04', 'P03'), 'VZ', to_rad(89.7080, unit='deg'))
    ]
    assert measurement_list == expected_measurements

def test_itersigma(non_standard_data):
    measurements = Measurements(non_standard_data, angle_unit='deg')
    sigma_list = list(measurements.itersigma())
    expected_sigma = [
        (('P01', 'P02'), 'sHD', 0.0024),
        (('P01', 'P02'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P01', 'P03'), 'sHD', 0.0027),
        (('P01', 'P03'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P02', 'P01'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P02', 'P03'), 'sHD', 0.0024),
        (('P02', 'P03'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P02', 'P04'), 'sHD', 0.0025),
        (('P02', 'P04'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P03', 'P01'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P03', 'P04'), 'sHD', 0.0023),
        (('P03', 'P04'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P04', 'P01'), 'sHD', 0.0026),
        (('P04', 'P01'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P04', 'P02'), 'sHZ', to_rad(0.0018, unit='deg')),
        (('P04', 'P03'), 'sHZ', to_rad(0.0018, unit='deg'))
    ]
    assert sigma_list == expected_sigma

def test_to_disp(non_standard_data):
    measurements = Measurements(non_standard_data, angle_unit='deg')
    converted = measurements.to_disp(angle_unit='deg')
    assert pytest.approx(converted['HZ'][0]) == 277.8982
    assert pytest.approx(converted['VZ'][0]) == 90.5143

def test_copy_with_type(non_standard_data):
    measurements = Measurements(non_standard_data)
    measurements_copy = measurements.copy_with_type()
    assert isinstance(measurements_copy, Measurements)
    assert id(measurements) != id(measurements_copy)
    pd.testing.assert_frame_equal(measurements, measurements_copy)
