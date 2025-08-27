# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.


import numpy as np
import pandas as pd

from ._constants import INVALID_INDEX
from .adjustment_results import AdjustmentResults


class Results(AdjustmentResults):
    def _rename_level(
        self, index: pd.MultiIndex, name: str, pos: int, inplace: bool = False
    ):
        names = list(index.names)
        names[pos] = name
        if not inplace:
            return index.set_names(names, inplace=False)
        index.set_names(names, inplace=True)

    def _get_obs_index(self) -> pd.MultiIndex:
        """Return MultiIndex object to describe observation adjustment results."""
        measurements = self._solver.dataset.measurements

        index = measurements.stack(future_stack=True).index
        return self._rename_level(index, "column", -1)

    def _get_matrix_coord_index(self) -> pd.MultiIndex:
        """Return MultiIndex object to describe coordinate indices in adjustment matrices."""
        coordinate_indices = self._solver.matrices.indexer.coordinate_indices

        mask = coordinate_indices != INVALID_INDEX
        filtered_indices = coordinate_indices[mask]

        index = filtered_indices.stack(future_stack=True).index
        return self._rename_level(index, "column", -1)

    def _get_matrix_orientation_index(self) -> pd.MultiIndex | None:
        """Return MultiIndex object to describe oreientation indices in adjustment matrices."""
        orientation_indices = self._solver.matrices.indexer.orientation_indices

        if orientation_indices is None:
            return

        mask = orientation_indices != INVALID_INDEX
        filtered_indices = orientation_indices[mask]

        return pd.MultiIndex.from_arrays(
            [filtered_indices.index, ["orientation"] * filtered_indices.size]
        )

    def _get_matrix_index(self) -> pd.MultiIndex:
        """Return MultiIndex object to describe adjustment matrices."""
        if self._matrix_orientation_index is not None:
            return self._matrix_coordinate_index.append(self._matrix_orientation_index)
        return self._matrix_coordinate_index

    def _get_adjusted_coordinate_sigmas(self) -> pd.DataFrame:
        """Return DataFrame containing adjusted coordinate sigmas."""
        coord_index = self._solver.dataset.controls.index
        coord_columns = self._solver.dataset.controls.coordinate_columns

        idx = self._solver.matrices.indexer.coordinate_indices
        mask = self._solver.matrices.indexer.coordinate_mask

        coord_var = self.covariance_X.values.diagonal()[idx[mask]]
        coord_sigmas = np.sqrt(coord_var)
        data = coord_sigmas.reshape(len(coord_index), len(coord_columns))

        return pd.DataFrame(
            data, index=coord_index, columns=[f"s{col}" for col in coord_columns]
        )

    def _get_coord_error_elipses(self) -> pd.DataFrame | None:
        """Return DataFrame containing adjusted coordinate error ellipses parameters."""
        coord_columns = self._solver.dataset.controls.coordinate_columns
        if not any(col in ["x", "y"] for col in coord_columns):
            return

        coord_idx = self._solver.dataset.controls.index
        cov_matrices = np.empty((len(coord_idx), 2, 2))

        for i, idx in enumerate(coord_idx):
            block_index = pd.MultiIndex.from_product([[idx], ["x", "y"]])
            cov_matrices[i] = self.covariance_X.loc[block_index, block_index].to_numpy()

        var_x = cov_matrices[:, 0, 0]
        var_y = cov_matrices[:, 1, 1]
        cov_xy = cov_matrices[:, 0, 1]

        term_1 = (var_x + var_y) / 2
        term_2 = np.sqrt(((var_x - var_y) / 2) ** 2 + cov_xy**2)

        a_sq = term_1 + term_2
        b_sq = term_1 - term_2
        phi = np.arctan2(2 * cov_xy, var_x - var_y) / 2

        return pd.DataFrame(
            {
                "a": np.sqrt(a_sq),
                "b": np.sqrt(b_sq),
                "phi": np.mod(phi, np.pi),
            },
            index=coord_idx,
        )

    def _get_adjusted_observation_values(self) -> pd.DataFrame:
        """Return DataFrame containing adjusted observation (measurement) values."""
        df = self.obs_residuals.rename(columns=lambda col: col[1:])
        return self._solver.dataset.measurements - df

    def _get_adjusted_observation_sigmas(self) -> pd.DataFrame:
        """Return DataFrame containing adjusted observation (measurement) sigmas."""
        data = self.covariance_Y.values.diagonal()
        df = pd.Series(data, index=self.covariance_Y.index).unstack(
            "column", sort=False
        )
        df.columns = [f"s{col}" for col in df.columns]
        return np.sqrt(df)
