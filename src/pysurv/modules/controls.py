# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pandas as pd


class Controls(pd.DataFrame):
    """
    A specialized DataFrame subclass for storing and manipulating control point data, including coordinates (x, y, z)
    and optionally their standard deviations (sx, sy, sz).

    --------------------------------------------------------------------------------------------------------------------
    Properties:

    - coord_label: (pd.Index): List of only coordinates labels.
    - labels: (pd.Index): List of coordinate labels and orientation constant, if in the set.
    - sigma_labels: (pd.Index): List of only coordinates sigma labels.
    - points_index: (dict): A dictionary mapping point labels (index) to their numerical positions.
    - coord_index: (dict): A dictionary mapping coordinate labels (index) to their numerical positions.
    - n_coords: (int): Number of control points coordinates and orietation constant if exist.

    --------------------------------------------------------------------------------------------------------------------
    Methods:

    - swap_xy(): Swap the x and y coordinates and their corresponding sigma values, if present in the DataFrame.
    - copy_with_type(): Create and return a deep copy of the Controls dataset.
    """
    acceptable_labels = {
        'coordinates': ['x', 'y', 'z'],
        'sigma': ['sx', 'sy', 'sz']
    }

    def __init__(self, data: pd.DataFrame, swap_xy: bool = False, *args, **kwargs):
        """
        Initialize the Controls DataFrame by standardizing column names and optionally swapping x and y coordinates.

        ----------------------------------------------------------------------------------------------------------------
        Arguments:
        - data: (pd.DataFrame): The input DataFrame containing control point data.
        - swap_xy: (bool): If True, swaps the x and y coordinates (and their corresponding sigma values if present).
        - *args, **kwargs: Additional positional and keyword arguments passed to the pandas DataFrame initializer.
        """
        # Standardize column names and set 'ID' as the index
        data = data.rename(columns={
            'id': 'ID', 'nr': 'ID',  # IDs
            'E': 'y', 'sE': 'sy',  # Easting to y
            'N': 'x', 'sN': 'sx',  # Northing to x
            'el': 'z', 'sel': 'sz',  # Elevation to z
            'h': 'z', 'sh': 'sz',  # Height to z
            'H': 'z', 'sH': 'sz'  # Height to z
        })
        # Convert point labels to string and set as index
        data['ID'] = data['ID'].astype(str)
        data.set_index(['ID'], inplace=True, verify_integrity=True)
        # Retain only acceptable columns for coordinates and sigma
        data = data.loc[:, data.columns.isin(
            Controls.acceptable_labels['coordinates'] +
            Controls.acceptable_labels['sigma']
        )]
        # Initialize the pandas DataFrame
        super().__init__(data, *args, **kwargs)

        if swap_xy:
            self.swap_xy()

    @property
    def coord_labels(self):
        """
        Returns: pd.Index: List of only coordinates labels
        """
        return self.columns[self.columns.isin(Controls.acceptable_labels['coordinates'])]

    @property
    def sigma_labels(self):
        """
        Returns: pd.Index: List of only coordinates sigma labels
        """
        return self.columns[self.columns.isin(Controls.acceptable_labels['sigma'])]

    @property
    def labels(self):
        """
        Returns: pd.Index: List of coordinates labels and orientation constant if exist
        """
        return self.columns[~self.columns.isin(Controls.acceptable_labels['sigma'])]

    @property
    def n_coords(self):
        """
        Returns: int: Number of control points coordinates and orietation constant if exist.
        """
        return self[self.labels].count().sum()

    @property
    def points_index(self):
        """
        Returns: dict: A dictionary mapping point labels (index) to their numerical positions.
        """
        return {point: idx for idx, point in enumerate(self.index)}

    @property
    def coord_index(self):
        """
        Returns: dict: A dictionary mapping coordinate labels (index) to their numerical positions.
        """
        controls_notna = self[self.labels].notna()
        coord_idx = dict()
        idx = 0
        for point, row in controls_notna.iterrows():
            coord_idx.update({point: {label: idx + i for i, label in
                                  enumerate(controls_notna.columns[controls_notna.loc[point]])}})
            idx += row.sum()
        return coord_idx

    def swap_xy(self):
        """
        Swap the x and y coordinates and their corresponding sigma values, if present in the DataFrame.
        """
        if {'x', 'y'}.issubset(self.columns):
            self[['x', 'y']] = self[['y', 'x']]
            if {'sx', 'sy'}.issubset(self.columns):
                self[['sx', 'sy']] = self[['sy', 'sx']]

    def copy_with_type(self):
        """
        Create and return a deep copy of the Controls dataset, preserving the Controls class type.

        ----------------------------------------------------------------------------------------------------------------
        Returns:

        - Controls: A copy of the current Controls object with its type preserved.
        """
        return Controls(self.copy().reset_index())
