# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
This module provides the `Project` class, which manages the setup and adjustment of surveying control network.

The `Project` class is responsible for:
- Importing controls and measurements data.
- Validating the adjustment methods provided.
- Creating an observation equations system.
- Performing the adjustment computation using the least squares method.
- Export report after computations.
"""
import os
import warnings
from pysurv import customizations
from pysurv.modules import Measurements
from pysurv.modules import importer
from pysurv.modules import adjustment
from pysurv.modules.basic import to_rad
from pysurv.modules.report import Report


class Project:
    """
    A class to represent a surveying project that manages control points, measurements, and adjustment computations.

    --------------------------------------------------------------------------------------------------------------------
    Attributes:

    - path: (str): The directory path where project files are stored.
    - name: (str): The name of the project, derived from the project path.
    - controls: (Controls): The control points dataset.
    - measurements: (Measurements): The measurements' dataset.
    - methods: (dict): The methods to be used for adjustment.
    - default_measurement_errors: (dict): Default measurement errors, with angle errors converted to radians.
    - iterate_params: (dict): Parameters for the iterative process during adjustment computations.
    - report_params: (dict): Parameters used for creating the adjustment results report.
    - observation_equations: (dict or None): The matrices representing system of observation equations.
      At the time of initialization, it takes the value None.
    - adjustment_results: (dict or None): The results of the adjustment process.
      At the time of initialization, it takes the value None.
    - report: (Report or None): The report object containing the results of the adjustment process.
      At the time of initialization, it takes the value None.

    --------------------------------------------------------------------------------------------------------------------
    Methods:

    - __str__(): Returns the name of the project.
    - create_observation_system(): Creates the system of observation equations.
    - adjust(): Performs the adjustment process on the dataset.
    """

    def __init__(self, path: str, methods: str = 'default', measurement_errors: str = 'default', swap_xy: bool = False,
                 angle_unit: str = 'grad', iterate_params: str = 'default', report_params: str = 'default'):
        """
        Initializes the Project with the provided path and configuration parameters.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - path: (str): The directory path to the project.
        - methods: (str, optional): Dictionary from customizations with the method set to use for adjustment.
          Defaults to 'default'.
        - measurement_errors: (str, optional): Dictionary from customizations with the default measurement sigma.
          Defaults to 'default'.
        - swap_xy: (bool, optional): Whether to swap the x and y coordinates when importing
          control points coordinates. Defaults to False.
        - angle_unit: (str, optional): The unit of the angles (e.g., 'grad', 'degree') in the measurements' dataset.
          Defaults to 'grad'.
        - iterate_params: (str, optional): Dictionary from customizations with parameters for the iterative process
          during adjustment computations. Defaults to 'default'.
        - report_params: (str, optional): Dictionary from customizations with parameters used for creating
          the adjustment results report. Defaults to 'default'.
        """
        # Validate and initialize the methods dictionary
        _check_methods(customizations.methods[methods])
        self.methods = customizations.methods[methods]
        # Convert angle errors to grad
        measurement_errors_rad = customizations.default_measurement_sigma[measurement_errors]
        for angle_sigma in Measurements.acceptable_labels['angular_measurements_sigma']:
            measurement_errors_rad[angle_sigma] = to_rad(measurement_errors_rad[angle_sigma], unit=angle_unit)
        self.default_measurement_errors = measurement_errors_rad
        # Initialize project attributes
        self.path = path
        self.name = os.path.basename(path)
        self.controls = importer.CSV.controls(os.path.join(path, 'controls.csv'), swap_xy=swap_xy)
        self.measurements = importer.CSV.measurements(os.path.join(path, 'measurements.csv'), angle_unit=angle_unit)
        self.iterate_params = customizations.iterate_params[iterate_params]
        self.report_params = customizations.report_params[report_params]
        self.observation_equations = None
        self.adjustment_results = None
        self.report = None

    def __str__(self):
        """
        Returns the name of the project.

        ----------------------------------------------------------------------------------------------------------------
        Returns:
        - str: The project name.
        """
        return f'project name: {self.name}'

    def create_observation_system(self):
        """
        Creates the system of observation equations.

        Observations equations are represented as matrices.
        This method updates the observation_equations attribute.

        ----------------------------------------------------------------------------------------------------------------
        Returns:
        - dict: The system of observation equations:
            - X: (np.ndarray): The coefficient matrix.
            - Y: (np.ndarray): The vector of observed values.
            - W: (np.ndarray, optional): The observations weight matrix.
            - R: (np.ndarray, optional): The matrix of network inner constraints.
            - sX: (np.ndarray, optional): The coefficient matrix for control points increments.
            - sW: (np.ndarray, optional): The weight matrix for the control points.
        """
        self.observation_equations = adjustment.matrices.equations_system(
            self.controls,
            self.measurements,
            methods=self.methods,
            default_sigma=self.default_measurement_errors
        )
        return self.observation_equations

    def adjust(self, report_format: str or None = None):
        """
        Performs the adjustment process on the dataset.

        If the observation equations have not been created, they are generated first.
        Adjustment computations are performed according to the selected methods.
        After the calculation is completed, it is possible to generate a report
        with the results to a file saved in the project folder.
        This method updates the adjustment_results and report attributes.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - report_format: (str, optional): The file format in which to export the report.
          If it is None, report will not be exported to the file. Defaults to None.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - dict: The results of the adjustment process:
            - n_iter: (int): The number of iterations performed.
            - sigma_zero: (list): List of sigma zero values from each iteration.
            - pt_sigma_zero: (list): List of control points increments sigma zero values from each iteration.
            - approx_coordinates: (Controls): The approximate control points coordinates before adjustment.
            - adjusted_coordinates: (Controls): The adjusted control points coordinates.
            - constraints_list: (list, optional): List of inner constraints applied to the network.
        """
        if self.observation_equations is None:
            self.create_observation_system()
        # Perform adjustment
        self.adjustment_results = adjustment.computations.adjust(
            self.controls,
            self.measurements,
            self.observation_equations,
            methods=self.methods,
            iterate_params=self.iterate_params
        )
        # Create report with results
        self.report = Report(self.adjustment_results,
                             self.measurements,
                             path=self.path,
                             methods=self.methods,
                             report_params=self.report_params
                             )
        # Export report to the particular format
        if report_format:
            report_format_function = self.report.formats.get(report_format)
            if report_format_function:
                report_format_function()
            else:
                warnings.warn(f'No report generated! Wrong file format in export_report parameter: {report_format}!')

        return self.adjustment_results


def _check_methods(methods: dict) -> bool:
    """
    Validates the provided adjustment methods against acceptable values.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - methods (dict): A dictionary containing the adjustment methods for 'observations' and 'free'.

    --------------------------------------------------------------------------------------------------------------------
    Returns:

    - bool: True if methods are valid, otherwise raises a ValueError.

    --------------------------------------------------------------------------------------------------------------------
    Raises:

    - ValueError: If the provided methods are not valid.
    """
    acceptable_methods = list(adjustment.robust.methods.keys()) + ['ordinary', 'weighted']

    if methods['observations'] not in acceptable_methods:
        raise ValueError(
            f"Invalid method for 'observations': {methods['observations']}. "
            f"Acceptable values are: {acceptable_methods}"
        )
    if methods['free'] is not None and methods['free'] not in acceptable_methods:
        raise ValueError(
            f"Invalid method for 'free': {methods['free']}. "
            f"Acceptable values are: None or {acceptable_methods}"
        )
    return True
