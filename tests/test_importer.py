# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest
from unittest.mock import patch, MagicMock
import numpy as np
import pandas as pd
from pysurv.modules.importer import CSV
from pysurv.modules import Controls, Measurements


@pytest.fixture
def controls_data():
    return pd.DataFrame(
                {'ID': ['P01', 'P02', 'P03', 'P04'],
                 'x': [317.08, 706.43, 1063.29, 852.32],
                 'y': [1606.89, 1785.03, 1552.50, 1293.80],
                 'z': [103.57, 99.73, 99.19, 97.48]
                 })


@pytest.fixture
def measurements_data():
    return pd.DataFrame(
                {'FROM': ['P01', 'P01', 'P02', 'P02', 'P02', 'P03', 'P03', 'P04', 'P04', 'P04'],
                 'TO': ['P02', 'P03', 'P01', 'P03', 'P04', 'P01', 'P04', 'P01', 'P02', 'P03'],
                 'HD': [428.173, 748.195, np.nan, 425.933, 512.430, np.nan, 333.814, 620.088, np.nan, np.nan],
                 'sHD': [0.0024, 0.0027, np.nan, 0.0024, 0.0025, np.nan, 0.0023, 0.0026, np.nan, np.nan],
                 'HZ': [308.7758, 276.8270, 384.1063, 120.0253, 75.1682, 27.2014, 88.2790, 270.0645, 222.1380, 160.2063],
                 'VZ': [100.5715, 100.3728, 99.4291, 100.0820, 100.2799, 99.6268, 100.3253, 99.3742, 99.7206, 99.6756]
                 })


@patch('pysurv.modules.importer.os.path.isfile')
@patch('pysurv.modules.importer.pd.read_csv')
def test_CSV_controls(mock_read_csv, mock_isfile, controls_data):
    mock_isfile.return_value = True
    mock_read_csv.return_value = controls_data

    with patch('pysurv.modules.importer.Controls') as mock_controls:
        mock_controls_instance = MagicMock(spec=Controls)
        mock_controls.return_value = mock_controls_instance

        result = CSV.controls('dummy_path.csv', swap_xy=True)

        mock_isfile.assert_called_once_with('dummy_path.csv')
        mock_read_csv.assert_called_once_with('dummy_path.csv')
        mock_controls.assert_called_once_with(controls_data, swap_xy=True)
        assert result == mock_controls_instance


@patch('pysurv.modules.importer.os.path.isfile')
@patch('pysurv.modules.importer.pd.read_csv')
def test_additional_args_and_kwargs_in_controls(mock_read_csv, mock_isfile, controls_data):
    mock_isfile.return_value = True
    mock_read_csv.return_value = controls_data

    with patch('pysurv.modules.importer.Controls.__init__') as mock_controls_init:
        mock_controls_init.return_value = None

        CSV.controls('dummy_path.csv', swap_xy=True, extra_arg1='value1', extra_kwarg='value2')

        mock_controls_init.assert_called_once_with(
            controls_data, swap_xy=True, extra_arg1='value1', extra_kwarg='value2'
        )


@patch('pysurv.modules.importer.os.path.isfile')
@patch('pysurv.modules.importer.pd.read_csv')
def test_measurements_import(mock_read_csv, mock_isfile, measurements_data):
    mock_isfile.return_value = True
    mock_read_csv.return_value = measurements_data

    with patch('pysurv.modules.importer.Measurements') as mock_measurements:
        mock_measurements_instance = MagicMock(spec=Measurements)
        mock_measurements.return_value = mock_measurements_instance

        result = CSV.measurements('dummy_path.csv', angle_unit='deg')

        mock_isfile.assert_called_once_with('dummy_path.csv')
        mock_read_csv.assert_called_once_with('dummy_path.csv')
        mock_measurements.assert_called_once_with(measurements_data, angle_unit='deg')
        assert result == mock_measurements_instance


@patch('pysurv.modules.importer.os.path.isfile')
@patch('pysurv.modules.importer.pd.read_csv')
def test_additional_args_and_kwargs_in_measurements(mock_read_csv, mock_isfile, measurements_data):
    mock_isfile.return_value = True
    mock_read_csv.return_value = measurements_data

    with patch('pysurv.modules.importer.Measurements.__init__') as mock_measurements_init:
        mock_measurements_init.return_value = None

        CSV.measurements('dummy_path.csv', angle_unit='deg', extra_arg1='value1', extra_kwarg='value2')

        mock_measurements_init.assert_called_once_with(
            measurements_data, angle_unit='deg', extra_arg1='value1', extra_kwarg='value2'
        )


@patch('pysurv.modules.importer.os.path.isfile')
def test_controls_import_invalid_path(mock_isfile):
    mock_isfile.return_value = False

    with pytest.raises(ValueError):
        CSV.controls('invalid_path.csv')


@patch('pysurv.modules.importer.os.path.isfile')
def test_measurements_import_invalid_path(mock_isfile):
    mock_isfile.return_value = False

    with pytest.raises(ValueError):
        CSV.measurements('invalid_path.csv')
