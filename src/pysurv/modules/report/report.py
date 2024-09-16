# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
This module provides the `Report` class, which store information about adjustment results and generate report file.

The `Report` class enables:
- Storing information about adjustment computations.
- Display report as string.
- Export report to the HTML file.
"""
import warnings
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from pysurv import customizations
from pysurv.modules.adjustment import robust
from pysurv.modules import Measurements
from pysurv.modules.basic import from_rad


class Report:
    """
    A class to generate adjustment report for a surveying control network.

    Report contains general information about adjustment computations, sigma evolution plots,
    control points coordinates, measurements and standard deviations.
    Information can be exported to string or to HTML file.

    --------------------------------------------------------------------------------------------------------------------
    Attributes:

    - approx_prec: (int): Precision for the approximate coordinates.
    - adj_prec: (int): Precision for the adjusted coordinates.
    - angle_unit: (str): Unit for angles ('grad', 'gon', or 'deg').
    - angle_precision: (int): Precision for angle values.
    - path: (str): The path where the report will be saved.
    - N: (np.ndarray): Matrix of normal equations.
    - cov_b: (np.ndarray): Covariance matrix of the unknowns (adjusted coordinates).
    - cov_Y: (np.ndarray): Covariance matrix of the measurements.
    - cov_r: (np.ndarray): Covariance matrix of the residuals.
    - approx: (Controls): Approximate coordinates.
    - adjusted: (Controls): Adjusted coordinates.
    - obs_residuals: (np.ndarray): Residuals of the observations.
    - obs_residuals_normalized: (np.ndarray): Normalized residuals.
    - obs_sigma_zero: (List[float]): List of observation sigma zero values per iteration.
    - pt_sigma_zero: (List[float]): List of control points' sigma zero values per iteration.
    - constraints: ([List[str]], optional): Applied inner constraints during the adjustment.
    - n_iter: (int): Number of iterations used in the adjustment.
    - n_measurements: (int): Number of measurements in the dataset.
    - n_fixed_coords: (int): Number of fixed reference coordinates.
    - n_sigma_coords: (int): Number of movable reference coordinates.
    - n_unknowns: (int): Number of coordinates and orientation constants to adjust.
    - n_constraints: (int): Number of inner constraints included to the control network.

    --------------------------------------------------------------------------------------------------------------------
    Methods:

    - _prepare_general_info(): Prepares a DataFrame containing general information about the adjustment process.
    - _create_sigma_plot(): Creates a plot showing the evolution of sigma zero values during the adjustment iterations.
    - _prepare_controls_table(): Prepares a DataFrame containing information about the control points after the adjustment.
    - _prepare_measurements_table(): Prepares a DataFrame containing information about the measurements after the adjustment.
    - to_string(): Generates a string representation of the report, optionally including the sigma plot.
    - html(): Generates an HTML representation of the report and saves it to the specified path.
    """
    def __init__(self, results: dict, measurements: Measurements, path: str,
                 methods=customizations.methods['default'],
                 report_params=customizations.report_params['default']):
        """
        Initializes the Report object. Prepare necessary tables with information about adjustment computations,
        controls and measurements.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - results: (dict): The results from the adjustment process, including matrices, coordinates, residuals, etc.
        - measurements: (Measurements): The measurements dataset used in the adjustment process.
        - path: (str): The path where the report will be saved.
        - methods: (dict): Methods and tuning constants used for adjustment.
        - report_params: (dict): Parameters for customizing the report format.
        """
        # --- UNPACK REPORT PARAMS ---
        self.approx_prec = report_params.get('approx_precision')
        self.adj_prec = report_params.get('adjusted_precision')
        self.angle_unit = report_params.get('angle_unit')
        self.angle_precision = report_params.get('angle_precision')
        self.path = path
        # --- UNPACK ADJUSTMENT RESULTS ---
        self.N = results.get('N')
        # Covariance matrices
        self.cov_b = results.get('cov_b')
        self.cov_Y = results.get('cov_Y')
        self.cov_r = results.get('cov_r')
        # Coordinates
        self.approx = results.get('approx')
        self.adjusted = results.get('adjusted')
        # Coordinate increments
        self.increments = (self.adjusted[self.adjusted.coord_labels] - self.approx[self.approx.coord_labels]) * 1000
        # Residuals
        self.obs_residuals = results.get('obs_residuals')
        with np.errstate(divide='ignore', invalid='ignore'):
            self.obs_sigma = np.sqrt(self.cov_r.diagonal())
        self.obs_residuals_normalized = results.get('r_norm')
        # Adjustment sigma
        self.obs_sigma_zero = results.get('sigma_zero')
        self.pt_sigma_zero = results.get('pt_sigma_zero')
        # Applied inner constraints, if any
        self.constraints = results.get('constraints')
        # General information
        self.n_iter = results.get('n_iter')
        self.n_measurements = results.get('n_measurements')
        self.n_fixed_coords = results.get('n_fixed_coords')
        self.n_sigma_coords = results.get('n_sigma_coords')
        self.n_unknowns = results.get('n_unknowns')
        self.n_constraints = len(self.constraints) if self.constraints is not None else 0
        # Prepare report components
        self.general_info = self._prepare_general_info(methods)
        self.plot = self._create_sigma_plot()
        self.controls_table = self._prepare_controls_table()
        self.measurements_table = self._prepare_measurements_table(measurements)
        # Footer with credits and contact information
        self.footer = [
            "Calculations made with pysurv, by Michal Predki. \n",
            "Github repository: https://github.com/mpredki99/pysurv.\n",
            "You can contact me on LinkedIn for feedback or if you discover a bug.\n",
            " "
        ]

        # Supported export formats
        self.formats = {
            'html': self.html
        }

    def _prepare_general_info(self, methods: dict) -> pd.DataFrame:
        """
        Prepare DataFrame containing general information about adjustment.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - methods: (dict): The methods and tuning constants values to be used for adjustment.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - pd.DataFrame: DataFrame with general information about adjustment.
        """
        # Prepare general information about adjustment
        general_info = {
            'Number of iterations': self.n_iter,
            'Observations deviation - obs_sigma_zero': self.obs_sigma_zero[-1],
            'Controls deviation - pt_sigma_zero': self.pt_sigma_zero[-1],
            'Number of measurements': self.n_measurements,
            'Number of applied inner constraints': self.n_constraints if self.n_constraints != 0 else pd.NA,
            'Number of fixed reference coordinates': self.n_fixed_coords if self.n_fixed_coords != 0 else pd.NA,
            'Number of movable reference coordinates': self.n_sigma_coords if self.n_sigma_coords != 0 else pd.NA,
            'Network defect': len(self.N) - np.linalg.matrix_rank(self.N),
            'Number of unknowns': self.n_unknowns,
            'Degrees of freedom': self.n_measurements - self.n_unknowns + self.n_constraints,
            'Method for adjust observations': methods['observations']
        }
        # Add specific information for robust methods if used in observations adjustment
        if methods['observations'] in robust.methods:
            general_info['Tuning constant for observations reweight'] = {
                0: methods['obs_c'] if methods['obs_c'] is not None else
                robust.theoretical_c[methods['observations']]
            }
        # Add specific information about free adjustment if applied
        if methods['free'] == 'ordinary':
            general_info['Method for free adjustment'] = 'pseudoinverse'
        elif methods['free'] in robust.methods or methods['free'] == 'weighted':
            general_info['Method for free adjustment'] = 'inner constraints'
            general_info['Applied constraints'] = ', '.join(self.constraints)
            if methods['free'] in robust.methods:
                general_info['Method for control points reweight'] = methods['free']
                general_info['Tuning constant for control points reweight'] = {
                    0: methods['free_c'] if methods['free_c'] is not None else
                    robust.theoretical_c[methods['free']]
                }

        return pd.DataFrame(general_info, index=[0]).T.dropna()

    def _create_sigma_plot(self):
        """
        Create a plot showing the evolution of sigma zero values during the iteration process.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - matplotlib.figure.Figure: The generated plot.
        """
        # Create plot of adjustment sigmas
        x = np.arange(1, len(self.obs_sigma_zero) + 1)
        fig = plt.figure(figsize=(10, 5))
        plt.title('Plot of adjustment sigma', fontfamily='cambria', fontweight='bold', fontsize=15)
        plt.plot(x, self.obs_sigma_zero, 'o-r', label='observations - σ₀')
        plt.plot(x, self.pt_sigma_zero, 'o-b', label='controls - σᵦ')
        plt.legend(prop=font_manager.FontProperties(family='Cambria'), fontsize=10)
        plt.xlabel('number of iteration', fontfamily='cambria', fontsize=12)
        plt.ylabel('σ values', fontfamily='cambria', fontsize=12)
        plt.xticks(x, fontfamily='cambria', fontsize=10)
        plt.yticks(fontfamily='cambria', fontsize=10)
        plt.ylim(0, np.max([np.max(self.pt_sigma_zero) * 1.1, np.max(self.obs_sigma_zero) * 1.1, 1.2]))
        plt.grid()
        return fig

    def _prepare_controls_table(self):
        """
        Prepare a DataFrame containing information about control points after the adjustment.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - pd.DataFrame: DataFrame with control points' results.
        """
        # Prepare controls table
        controls_sigma = pd.DataFrame()
        # Get the row and column indices of the non-NaN values
        row_indices, col_indices = np.where(~self.adjusted[self.adjusted.labels].isna())
        # Apply the increments to the corresponding positions
        for row, col, sigma in zip(row_indices, col_indices, self.cov_b.diagonal()):
            if self.adjusted.labels[col] == 'o':
                continue
            controls_sigma.loc[self.adjusted.index[row], self.adjusted.labels[col]] = sigma
        with np.errstate(divide='ignore', invalid='ignore'):
            controls_sigma = np.sqrt(controls_sigma) * 1000
        # Normalize the coords increments
        increments_normalized = self.increments.div(controls_sigma).where(controls_sigma != 0)
        # Populate the table with approximate, increments, and adjusted coordinates
        controls_table = pd.DataFrame()
        for coord_label in self.adjusted.coord_labels:
            controls_table[f'{coord_label}_0_ [m]'] = self.approx[coord_label].round(self.approx_prec)
            controls_table[f'd{coord_label} [mm]'] = self.increments[coord_label].round(self.adj_prec - 3)
            controls_table[f'{coord_label} [m]'] = self.adjusted[coord_label].round(self.adj_prec)
        # Populate the table with standard deviations
        for coord_label in self.adjusted.coord_labels:
            controls_table[f'_s_{coord_label} [mm]'] = controls_sigma[coord_label].round(self.adj_prec - 2)
        # Populate the table with normalized increments
        for coord_label in self.adjusted.coord_labels:
            controls_table[f'd{coord_label}/_s_{coord_label}'] = increments_normalized[coord_label].round(2)

        return controls_table

    def _prepare_measurements_table(self, measurements):
        """
        Prepare a DataFrame containing information about measurements after adjustment.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - measurements (Measurements): The measurements dataset used in the adjustment.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - pd.DataFrame: DataFrame containing information about measurements.
        """
        # Helper function to apply the necessary transformations
        def apply_transformations(col, precision, multiplier: float = 1):
            """
            Apply necessary transformations to the measurement columns.

            ------------------------------------------------------------------------------------------------------------
            Arguments:

            - col (str): The column name.
            - precision (int): Precision for rounding.
            - multiplier (float, optional): Multiplier to apply to the column values. Default is 1.
            """
            measurements_table[col] *= multiplier
            measurements_table[col] = measurements_table[col].astype(float).round(precision)

        measurements_table = pd.DataFrame()
        # Precompute empty fields
        for meas_type in measurements.types:
            measurements_table[f'{meas_type}_0_'] = measurements[meas_type]
            measurements_table[f'v{meas_type}'] = pd.Series(index=measurements.index)
            measurements_table[meas_type] = pd.Series(index=measurements.index)
            measurements_table[f'_s_{meas_type}'] = pd.Series(index=measurements.index)
            measurements_table[f'v/_s_v {meas_type}'] = pd.Series(index=measurements.index)
        # Assign values
        for i_meas, (index, meas_type, value) in enumerate(measurements.itermeasurements()):
            measurements_table.at[index, f'v{meas_type}'] = self.obs_residuals[i_meas]
            measurements_table.at[index, meas_type] = measurements.at[index, meas_type] + self.obs_residuals[i_meas]
            measurements_table.at[index, f'_s_{meas_type}'] = self.obs_sigma[i_meas]
            measurements_table.at[index, f'v/_s_v {meas_type}'] = self.obs_residuals_normalized[i_meas]
        # Apply transformations based on column type
        for column in measurements_table.columns:
            if column.startswith('v/_s_v '):
                apply_transformations(
                    column,
                    2
                )
            elif column in Measurements.angle_column_types:
                if column.startswith('v') or column.startswith('_s_'):
                    measurements_table[column] = from_rad(measurements_table[column], unit=self.angle_unit)
                    apply_transformations(
                        column,
                        self.angle_precision - (4 if column.startswith('v') else 3),
                        10000 if self.angle_unit == 'grad' else 3600
                    )
                else:
                    measurements_table[column] = from_rad(measurements_table[column], unit=self.angle_unit)
                    apply_transformations(column, self.angle_precision)
            else:
                if column.startswith('v') or column.startswith('_s_'):
                    apply_transformations(
                        column,
                        self.adj_prec - (3 if column.startswith('v') else 2),
                        1000
                    )
                elif column.endswith('_0_'):
                    apply_transformations(
                        column,
                        self.approx_prec
                    )
                else:
                    apply_transformations(
                        column,
                        self.adj_prec
                    )

        return measurements_table

    def to_string(self, show_plot: bool = True) -> str:
        """
        Generate a string representation of the report.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - show_plot: (bool, optional): If True, include the sigma plot in the string. Default is True.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - str: String representation of the report.
        """
        from ._to_string import to_string
        return to_string(self.plot,
                         self.general_info,
                         self.controls_table,
                         self.measurements_table,
                         self.footer,
                         show_plot=show_plot
                         )

    def html(self, path=None):
        """
        Generate an HTML representation of the report and export it to the file.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - path: (str, optional): The path to save the HTML file. If not provided,
          default is the path provided at initialization.
        """
        import os
        from ._to_html import to_html

        if path is None:
            path = self.path
        # Prepare content of the file
        html_content = to_html(self.plot,
                               self.general_info,
                               self.controls_table,
                               self.measurements_table,
                               self.footer
                              )
        # Write the HTML content to a file
        with open(os.path.join(path, 'report.html'), 'w') as file:
            file.write(html_content)
        print("HTML file created successfully!")
