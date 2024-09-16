# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest
import numpy as np
import pandas as pd
import os
from matplotlib.figure import Figure
from unittest.mock import patch, mock_open
from pysurv.modules import Measurements, Controls, Report

@pytest.fixture
def create_report():
    results = {
        'cov_b': np.eye(12),
        'cov_Y': np.eye(26),
        'cov_r': np.eye(26),
        'approx': Controls(
            pd.DataFrame(
                {'ID': ['P01', 'P02', 'P03', 'P04'],
                 'x': [317.08, 706.43, 1063.29, 852.32],
                 'y': [1606.89, 1785.03, 1552.50, 1293.80],
                 'z': [103.57, 99.73, 99.19, 97.48]
                 })),
        'adjusted': Controls(
            pd.DataFrame(
                {'ID': ['P01', 'P02', 'P03', 'P04'],
                 'x': [317.081, 706.431, 1063.291, 852.321],
                 'y': [1606.891, 1785.031, 1552.501, 1293.801],
                 'z': [103.571, 99.731, 99.191, 97.481]
                 })),
        'obs_residuals': np.array([-0.01, 0.02, -0.03, 0.04, -0.05, 0.06, -0.07, 0.08, -0.09, 0.10, -0.11, 0.12, -0.13, 0.14, -0.15, 0.16, -0.17, -0.18, -0.19, 0.20, -0.21, 0.22, -0.23, 0.24, -0.25, 0.26]),
        'r_norm': np.array([0.01, -0.02, 0.03, -0.04, 0.05, -0.06, 0.07, -0.08, 0.09, -0.10, 0.11, -0.12, 0.13, -0.14, 0.15, -0.16, 0.17, -0.18, 0.19, -0.20, 0.21, -0.22, 0.23, -0.24, 0.25, -0.26]),
        'sigma_zero': [1.5, 1.2, 1.0],
        'pt_sigma_zero': [0.5, 0.7, 1.0],
        'constraints': ['translation x', 'translation y', 'translation z'],
        'n_iter': 3,
    }
    measurements_data = pd.DataFrame(
                {'FROM': ['P01', 'P01', 'P02', 'P02', 'P02', 'P03', 'P03', 'P04', 'P04', 'P04'],
                 'TO': ['P02', 'P03', 'P01', 'P03', 'P04', 'P01', 'P04', 'P01', 'P02', 'P03'],
                 'HD': [428.173, 748.195, np.nan, 425.933, 512.430, np.nan, 333.814, 620.088, np.nan, np.nan],
                 'sHD': [0.0024, 0.0027, np.nan, 0.0024, 0.0025, np.nan, 0.0023, 0.0026, np.nan, np.nan],
                 'HZ': [308.7758, 276.8270, 384.1063, 120.0253, 75.1682, 27.2014, 88.2790, 270.0645, 222.1380, 160.2063],
                 'VZ': [100.5715, 100.3728, 99.4291, 100.0820, 100.2799, 99.6268, 100.3253, 99.3742, 99.7206, 99.6756]
                 }
        )
    measurements = Measurements(measurements_data)
    path = 'test'
    methods = {'observations': 'huber',
               'obs_c': 1.345,
               'free': 'tukey',
               'free_c': None,
               }
    report = Report(
        results=results,
        measurements=measurements,
        path=path,
        methods=methods,
    )
    return report

def test_general_info(create_report):
    report = create_report
    # Check if the actual shape matches the expected shape
    assert report.general_info.shape == (12, 1)
    # Test the general information is DataFrame
    assert isinstance(report.general_info, pd.DataFrame)

def test_sigma_plot(create_report):
    report = create_report
    # Test the sigma plot creation
    assert isinstance(report.plot, Figure)

def test_controls_table(create_report):
    report = create_report
    # Check if the actual shape matches the expected shape
    assert report.controls_table.shape == (4, 15)
    # Test the controls table is DataFrame
    assert isinstance(report.controls_table, pd.DataFrame)

def test_measurements_table(create_report):
    report = create_report
    # Check if the actual shape matches the expected shape
    assert report.measurements_table.shape == (10, 15)
    # Test the measurements table is DataFrame
    assert isinstance(report.measurements_table, pd.DataFrame)

def test_to_string(create_report):
    report = create_report
    # Test the string representation of the report
    report_string = report.to_string(show_plot=False)
    # Test the general information DataFrame
    assert isinstance(report_string, str)
    assert 'Number of iterations' in report_string
    assert 'Calculations made with pysurv' in report_string


@patch('pysurv.modules.report._to_html.to_html')
@patch('builtins.open', new_callable=mock_open)
def test_html_method(mock_file, mock_to_html, create_report, capsys):
    """
    Test the html method to verify that it creates a file in the right place with the right content,
    and that it informs the user that the report creation process is complete.
    """
    # Create instance of Report class
    report = create_report

    # Call the html method
    report.html()
    # Capture the printed information about export
    out, err = capsys.readouterr()

    # Verify that to_html was called correctly with the expected arguments
    mock_to_html.assert_called_once_with(
        report.plot,
        report.general_info,
        report.controls_table,
        report.measurements_table,
        report.footer
    )
    # Verify that the file was opened in write mode
    mock_file.assert_called_once_with(os.path.join(report.path, 'report.html'), 'w')
    # Verify the writer received the proper file content
    mock_file.return_value.write.assert_called_once_with(mock_to_html.return_value)
    # Verify the prompt was displayed correctly
    assert out.strip() == 'HTML file created successfully!'