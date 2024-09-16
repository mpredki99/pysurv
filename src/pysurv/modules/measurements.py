# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from typing import Tuple
import pandas as pd
import numpy as np
from pysurv.modules.basic import to_rad, from_rad


class Measurements(pd.DataFrame):
    """
    A specialized pd.DataFrame subclass for storing and manipulating measurements data,
    and optionally their standard deviations.

    --------------------------------------------------------------------------------------------------------------------
    Properties:

    - types: (pd.Index): List of only measurements types labels.
    - sigma: (pd.Index): List of only measurements sigma labels.
    - linear: (pd.Index): List of only linear measurements types labels.
    - linear_sigma: (pd.Index): List of only linear measurements sigma labels.
    - angular: (pd.Index): List of only angular measurements types labels.
    - angular_sigma: (pd.Index): List of only angular measurements sigma labels.
    - n_measurements: (int): Number of measurements.

    --------------------------------------------------------------------------------------------------------------------
    Methods:

    - itermeasurements(): Iterates through all the measurements in the dataset.
    - itersigma(): Iterates through all the sigmas in the dataset.
    - iteritems(): Iterates through all the measurements and sigma values in the dataset.
    - to_disp(): Returns a deep copy of the dataset with angles converted to the specified unit.
    - copy_with_type(): Create and return a deep copy of the Measurements dataset.
    """
    acceptable_labels = {
        'linear_measurements': ['SD', 'HD', 'VD', 'dx', 'dy', 'dz'],
        'linear_measurements_sigma': ['sSD', 'sHD', 'sVD', 'sdx', 'sdy', 'sdz'],
        'angular_measurements': ['A', 'HZ', 'VZ', 'VH'],
        'angular_measurements_sigma': ['sA', 'sHZ', 'sVZ', 'sVH'],
    }
    angle_column_types = [
        'A_0_', 'vA', 'A', '_s_A', 'sA',
        'HZ_0_', 'vHZ', 'HZ', '_s_HZ', 'sHZ',
        'VZ_0_', 'vVZ', 'VZ', '_s_VZ', 'sVZ',
        'VH_0_', 'vVH', 'VH', '_s_VH', 'sVH'
    ]

    def __init__(self, data: pd.DataFrame, angle_unit: str = 'grad', *args, **kwargs):
        """
        Initialize the Controls DataFrame by standardizing column names and convert angles to radians.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - data: (pd.DataFrame): The input DataFrame containing measurements data.
        - angle_unit: (str): Unit of angular measurements in the dataset ('grad', 'gon' or 'deg').
        - *args, **kwargs: Additional positional and keyword arguments passed to the pandas DataFrame initializer.

        ----------------------------------------------------------------------------------------------------------------
        Raises:

        - ValueError: If angle_unit is not 'grad', 'gon', 'deg' or 'rad'.
        """
        if angle_unit not in ['gon', 'grad', 'deg', 'rad']:
            raise ValueError('Invalid unit. Use "grad", "gon", "deg" or "rad".')
        # Standardize column names and set ('FROM', 'TO') as the MultiIndex
        data = data.rename(columns={
            'from': 'FROM', 'to': 'TO',
            'sd': 'SD', 'hd': 'HD', 'vd': 'VD', 'DX': 'dx', 'DY': 'dy', 'DZ': 'dz',
            'a': 'A', 'hz': 'HZ', 'vz': 'VZ', 'vh': 'VH',
            'ssd': 'sSD', 'shd': 'sHD', 'svd': 'sVD', 'sDX': 'sdx', 'sDY': 'sdy', 'sDZ': 'sdz',
            'sa': 'sA', 'shz': 'sHZ', 'svz': 'sVZ', 'svh': 'sVH',
        })
        # Convert point labels to string and set as index
        data[['FROM', 'TO']] = data[['FROM', 'TO']].astype(str)
        data.set_index(['FROM', 'TO'], inplace=True, verify_integrity=True)
        # Retain only acceptable columns for measurements and sigma
        data = data.loc[:, data.columns.isin(Measurements.acceptable_labels['linear_measurements'] +
                                             Measurements.acceptable_labels['linear_measurements_sigma'] +
                                             Measurements.acceptable_labels['angular_measurements'] +
                                             Measurements.acceptable_labels['angular_measurements_sigma'])]
        # Convert angles to radians
        if angle_unit != 'rad':
            angle_columns = (Measurements.acceptable_labels['angular_measurements'] +
                             Measurements.acceptable_labels['angular_measurements_sigma'])
            angle_columns = [col for col in angle_columns if col in data.columns]
            if angle_columns:
                data[angle_columns] = to_rad(data[angle_columns], unit=angle_unit)
        # Initialize the pands DataFrame
        super().__init__(data, *args, **kwargs)

    @property
    def types(self):
        """
        Returns: pd.Index: List of only measurements types labels.
        """
        return self.columns[self.columns.isin(
            Measurements.acceptable_labels['linear_measurements'] +
            Measurements.acceptable_labels['angular_measurements'])
        ]

    @property
    def sigma(self):
        """
        Returns: pd.Index: List of only measurements sigma labels.
        """
        return self.columns[self.columns.isin(
            Measurements.acceptable_labels['linear_measurements_sigma'] +
            Measurements.acceptable_labels['angular_measurements_sigma'])
        ]

    @property
    def linear(self):
        """
        Returns: pd.Index: List of only linear measurements types labels.
        """
        return self.columns[self.columns.isin(Measurements.acceptable_labels['linear_measurements'])]

    @property
    def linear_sigma(self):
        """
        Returns: pd.Index: List of only linear measurements sigma labels.
        """
        return self.columns[self.columns.isin(Measurements.acceptable_labels['linear_measurements_sigma'])]

    @property
    def angular(self):
        """
        Returns: pd.Index: List of only angular measurements types labels.
        """
        return self.columns[self.columns.isin(Measurements.acceptable_labels['angular_measurements'])]

    @property
    def angular_sigma(self):
        """
        Returns: pd.Index: List of only angular measurements sigma labels.
        """
        return self.columns[self.columns.isin(Measurements.acceptable_labels['angular_measurements_sigma'])]

    @property
    def n_measurements(self):
        """
        Returns: int: Number of measurements.
        """
        return self[self.types].count().sum()

    def itermeasurements(self, show_empty: bool = False) -> Tuple[pd.MultiIndex, str, float]:
        """
        Iterates through all the measurements in the dataset.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - show_empty: (bool, optional): If True, yields measurements with NaN values. If False, skips NaN values.
          Default is False.

        ----------------------------------------------------------------------------------------------------------------
        Yields:

        - tuple: Measurements in dataset:
            * (pd.MultiIndex): Identifier of points 'FROM' and 'TO'.
            * (str): The measurement type.
            * (float): The measurement value.
        """
        for index, row in self.iterrows():
            for measurement_type in self.types:
                value = row.get(measurement_type)
                if np.isnan(value) and not show_empty:
                    continue
                yield index, measurement_type, value

    def itersigma(self, show_empty: bool = False) -> Tuple[pd.MultiIndex, str, float]:
        """
        Iterates through all the sigmas in the dataset.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - show_empty: (bool, optional): If True, yields measurements with NaN values. If False, skips NaN values.
          Default is False.

        ----------------------------------------------------------------------------------------------------------------
        Yields:

        - tuple: A tuple containing the index tuple:
            * (pd.MultiIndex): Identifier of points 'FROM' and 'TO'.
            * (str): The sigma type.
            * (float): The sigma value.
        """
        for index, row in self.iterrows():
            for sigma_type in self.sigma:
                value = row.get(sigma_type)
                if np.isnan(value) and not show_empty:
                    continue
                yield index, sigma_type, value

    def iteritems(self, show_empty: bool = False) -> Tuple[pd.MultiIndex, str, float]:
        """
        Iterates through all the measurements and sigma values in the dataset.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - show_empty: (bool, optional): If True, yields measurements with NaN values. If False, skips NaN values.
          Default is False.

        ----------------------------------------------------------------------------------------------------------------
        Yields:

        - tuple: A tuple containing the index tuple
            * (pd.MultiIndex): Identifier of points 'FROM' and 'TO'.
            * (str): The column type.
            * (float): The value.
        """
        for index, row in self.iterrows():
            for item_type in self.columns:
                value = row.get(item_type)
                if np.isnan(value) and not show_empty:
                    continue
                yield index, item_type, value

    def to_disp(self, angle_unit: str = 'grad'):
        """
        Returns a deep copy of the dataset with angles converted to the specified unit.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:

        - angle_unit: (str, optional): The unit in which angles should be displayed.
          Supported units include 'grad', 'gon' and 'deg'. Default is 'grad'.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - Measurements: A copy of the current dataset with angles converted to the specified unit.
        """
        dataset = self.copy_with_type()
        angle_columns = [col for col in self.columns if col in Measurements.angle_column_types]
        if angle_columns:
            dataset[angle_columns] = from_rad(self[angle_columns], angle_unit)

        return dataset

    def copy_with_type(self):
        """
        Create and return a deep copy of the Measurements DataFrame, preserving the Measurements class type.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - Measurements: A copy of the current Measurements object with its type preserved.
        """
        return Measurements(self.copy().reset_index(), angle_unit='rad')
