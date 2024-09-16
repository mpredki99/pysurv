# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest
import numpy as np
import pandas as pd
from pysurv.modules import Controls

@pytest.fixture
def standard_data():
    # Sample data for testing
    return pd.DataFrame(
            {'ID': ['P01', 'P02', 'P03', 'P04'],
             'x': [317.08, 706.43, 1063.29, 852.32],
             'y': [1606.89, 1785.03, 1552.50, 1293.80],
             'z': [103.57, 99.73, 99.19, 97.48],
             'sx': [0.01, 0.01, 0.01, np.nan],
             'sy': [0.01, 0.01, 0.01, np.nan],
             'sz': [0.05, 0.05, 0.05, np.nan]
             })

@pytest.fixture
def non_standard_data():
    # Sample data for testing
    return pd.DataFrame(
            {'id': ['P01', 'P02', 'P03', 'P04'],
             'N': [317.08, 706.43, 1063.29, 852.32],
             'E': [1606.89, 1785.03, 1552.50, 1293.80],
             'H': [103.57, 99.73, 99.19, 97.48],
             'sN': [0.01, 0.01, 0.01, np.nan],
             'sE': [0.01, 0.01, 0.01, np.nan],
             'sH': [0.05, 0.05, 0.05, np.nan]
             })

def test_initialization_standard_columns(standard_data):
    controls = Controls(standard_data)
    # Check that dataset is required type
    assert isinstance(controls, Controls)
    # Check that columns match expected labels
    pd.testing.assert_index_equal(controls.columns, pd.Index(['x', 'y', 'z', 'sx', 'sy', 'sz']))
    # Ensure the data is properly initialized
    pd.testing.assert_frame_equal(controls, standard_data.set_index('ID'))

def test_initialization_with_nonstandard_columns(standard_data, non_standard_data):
    controls = Controls(non_standard_data)
    # Check that dataset is required type
    assert isinstance(controls, Controls)
    # Verify that columns were renamed properly
    pd.testing.assert_index_equal(controls.columns, pd.Index(['x', 'y', 'z', 'sx', 'sy', 'sz']))
    # Ensure the data is properly initialized
    pd.testing.assert_frame_equal(controls, standard_data.set_index('ID'))

def test_swap_xy(standard_data, non_standard_data):
    controls = Controls(non_standard_data, swap_xy=True)
    # Check that x and y coordinates and their sigmas were swapped
    swapped_data = standard_data.set_index('ID')
    pd.testing.assert_series_equal(controls['x'].rename('y'), swapped_data['y'])
    pd.testing.assert_series_equal(controls['y'].rename('x'), swapped_data['x'])
    pd.testing.assert_series_equal(controls['sx'].rename('sy'), swapped_data['sy'])
    pd.testing.assert_series_equal(controls['sy'].rename('sx'), swapped_data['sx'])

def test_copy_with_type(standard_data):
    controls = Controls(standard_data)
    controls_copy = controls.copy_with_type()
    # Check that the copied object is indeed a Controls instance
    assert isinstance(controls_copy, Controls)
    # Check that the data is identical but not the same object (deep copy)
    assert id(controls) != id(controls_copy)
    pd.testing.assert_frame_equal(controls, controls_copy)

def test_coord_labels_property(non_standard_data):
    controls = Controls(non_standard_data)
    assert list(controls.coord_labels) == ['x', 'y', 'z']

def test_sigma_labels_property(non_standard_data):
    controls = Controls(non_standard_data)
    assert list(controls.sigma_labels) == ['sx', 'sy', 'sz']

def test_labels_property(non_standard_data):
    controls = Controls(non_standard_data)
    assert list(controls.labels) == ['x', 'y', 'z']

def test_map_index(non_standard_data):
    controls = Controls(non_standard_data)
    expected_map = {'P01': 0, 'P02': 1, 'P03': 2, 'P04': 3}
    assert controls.map_index == expected_map

def test_invalid_column_names():
    invalid_data = pd.DataFrame(
        {'id': ['P01', 'P02', 'P03', 'P04'],
         'axis_n': [317.08, 706.43, 1063.29, 852.32],
         'axis_e': [1606.89, 1785.03, 1552.50, 1293.80],
         'time': [103.57, 99.73, 99.19, 97.48],
         'std_n': [0.01, 0.01, 0.01, np.nan],
         'std_e': [0.01, 0.01, 0.01, np.nan],
         'std_t': [0.05, 0.05, 0.05, np.nan]
         })
    controls = Controls(invalid_data)
    # Check that no columns are retained as none are valid
    assert controls.empty
