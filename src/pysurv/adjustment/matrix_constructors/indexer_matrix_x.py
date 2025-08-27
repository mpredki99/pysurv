# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from functools import cached_property

import numpy as np
import pandas as pd

from pysurv.data.dataset import Dataset

from .._constants import INVALID_INDEX


class IndexerMatrixX:
    """
    Indexer for mapping control point coordinates and station orientations to matrix X column indices.

    This class provides methods to lazy generation and access index mappings for control point coordinates
    and station orientations, which are used as columns in the design matrix (X)for least squares
    adjustment computations.
    """

    def __init__(self, dataset: Dataset) -> None:
        self._controls = dataset.controls
        self._stations = dataset.stations

    @cached_property
    def coordinate_indices(self) -> pd.DataFrame:
        """Return a DataFrame mapping control points to their columns in matrix X."""
        return self._map_coordinate_index()

    @cached_property
    def coordinate_mask(self) -> pd.DataFrame:
        """Return a DataFrame mask that shows if coordinate index is valid value."""
        return self.coordinate_indices != INVALID_INDEX

    @cached_property
    def filtered_coordinate_indices(self) -> pd.DataFrame:
        """Return a coordinate_indices slice that contains only valid indices."""
        return self.coordinate_indices[self.coordinate_mask]

    @cached_property
    def orientation_indices(self) -> pd.Series:
        """Return a Series mapping station orientations to their columns in matrix X."""
        if "orientation" not in self._stations.columns:
            return

        if self.coordinate_indices is not None:
            return self._map_orientation_index()

    @cached_property
    def orientation_mask(self) -> pd.DataFrame:
        """Return a DataFrame mask that shows if orienatation index is valid value."""
        if self.orientation_indices is None:
            return
        return self.orientation_indices != INVALID_INDEX

    @cached_property
    def filtered_orientation_indices(self) -> pd.DataFrame:
        """Return a orientation_indices slice that contains only valid indices."""
        if self.orientation_indices is None:
            return
        return self.orientation_indices[self.orientation_mask]

    def _map_coordinate_index(self) -> None:
        """Map control point coordinates to their corresponding column indices in matrix X."""
        coord_idx = self._controls.index
        coord_columns = self._controls.coordinate_columns

        coord_values = self._controls.coordinates.values
        mask = np.isfinite(coord_values)

        index_map = np.full(coord_values.shape, INVALID_INDEX, dtype=int)
        index_map[mask] = np.arange(np.count_nonzero(mask))

        return pd.DataFrame(index_map, index=coord_idx, columns=coord_columns)

    def _map_orientation_index(self) -> None:
        """Map station orientations to their corresponding column indices in matrix X."""
        orientations = self._stations.orientation.dropna()

        start_idx = self.coordinate_indices.max().max() + 1
        end_idx = start_idx + len(orientations)

        index_values = np.arange(start_idx, end_idx)

        return pd.Series(index_values, index=orientations.index, name="orientation_idx")
