# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
This module provides functionalities for import controls and measurements datasets.

Imported datasets are used to create instance of Controls and Measurements classes.
"""
import os
import pandas as pd
from .controls import Controls
from .measurements import Measurements


class CSV:
    """
    Class used to import datasets from CSV files.

    --------------------------------------------------------------------------------------------------------------------
    Methods:

    - controls(): Imports a CSV file containing control points coordinates and sigma values.
    - measurements(): Imports a CSV file containing measurements and sigma values.
    """
    @staticmethod
    def controls(path: str, swap_xy: bool = False, *args, **kwargs) -> Controls:
        """
        Imports a CSV file containing control points coordinates and sigma values.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - path: (str): Path to the CSV file containing the controls' dataset.
        - swap_xy: (bool): Whether to swap the values of x and y coordinates. Defaults to False.
        - *args, **kwargs: Additional positional and keyword arguments passed to the pandas DataFrame initializer.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - Controls: An instance of the Controls class containing the controls dataset.

        ----------------------------------------------------------------------------------------------------------------
        Raises:

        - ValueError: If the provided path does not point to a valid file.
        """
        if not os.path.isfile(path):
            raise ValueError(f'Invalid path to the controls file: {path}')

        # Read CSV file and create Controls instance
        data = pd.read_csv(path)
        return Controls(data, swap_xy=swap_xy, *args, **kwargs)

    @staticmethod
    def measurements(path: str, angle_unit: str = 'grad', *args, **kwargs) -> Measurements:
        """
        Imports a CSV file containing measurements and sigma values.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - path: (str): Path to the CSV file containing the measurements' dataset.
        - angle_unit: (str): Unit of angular measurements in the dataset. Defaults to 'grad'.
        - *args, **kwargs: Additional positional and keyword arguments passed to the pandas DataFrame initializer.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - Measurements: An instance of the Measurements class containing the measurements dataset.

        ----------------------------------------------------------------------------------------------------------------
        Raises:

        - ValueError: If the provided path does not point to a valid file.
        """
        if not os.path.isfile(path):
            raise ValueError(f'Invalid path to the measurements file: {path}')

        # Read CSV file and create Measurements instance
        data = pd.read_csv(path)
        return Measurements(data, angle_unit=angle_unit, *args, **kwargs)
