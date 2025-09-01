# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from abc import ABC, abstractmethod

import pandas as pd

from pysurv.data.controls import Controls
from pysurv.utils.utils import refreshable_property


class AdjustmentResults(ABC):
    def __init__(self, solver):
        self._solver = solver

        self._observation_index = self._get_obs_index()
        self._matrix_coordinate_index = self._get_matrix_coord_index()
        self._matrix_orientation_index = self._get_matrix_orientation_index()
        self._matrix_index = self._get_matrix_index()
        
    def __bool__(self):
        return self.n_iter > 0

    @property
    def solver(self):
        return self._solver

    @property
    def _is_deprecated(self):
        return self._solver.current_iter > self.__dict__.get("n_iter", 0)

    @property
    def dataset(self):
        return self._solver.dataset

    @property
    def approx_coordinates(self) -> Controls:
        return self._solver.approx_coordinates.rename(
            columns={col: f"{col}_0" for col in self._solver.approx_coordinates.columns}
        )

    @property
    def inner_constraints(self) -> list:
        if self._solver.methods.free_adjustment == "ordinary":
            return ["pseudoinverse"]
        return self._solver.matrices.inner_constraints

    @property
    def methods(self):
        return self._solver.methods

    @property
    def n_measurements(self):
        return self._solver.matrices.n_measuremetns

    @property
    def n_unknowns(self):
        return self._solver.matrices.n_unknowns

    @property
    def degrees_of_freedom(self):
        return self._solver.matrices.degrees_of_freedom

    @property
    def n_coord_corrections(self):
        return self._solver.n_coord_corrections

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def n_iter(self) -> int:
        return self._solver.current_iter

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def matrix_G(self) -> pd.DataFrame | None:
        """Return G matrix as indexed DataFrame."""
        if self._solver.matrix_G is None:
            return

        return pd.DataFrame(
            self._solver.matrix_G,
            index=self._matrix_index,
            columns=self._matrix_index,
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def inv_matrix_G(self) -> pd.DataFrame | None:
        """Return inverse of G matrix as indexed DataFrame."""
        if self._solver.inv_matrix_G is None:
            return

        return pd.DataFrame(
            self._solver.inv_matrix_G,
            index=self._matrix_index,
            columns=self._matrix_index,
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def cross_product(self) -> pd.DataFrame | None:
        """Rerurn cross product"""
        if self._solver.cross_product is None:
            return

        return pd.Series(
            self._solver.cross_product,
            index=self._matrix_index,
            name="cross_product",
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def covariance_X(self) -> pd.DataFrame | None:
        """Return the covariance matrix of X."""
        if self._solver.covariance_X is None:
            return

        return pd.DataFrame(
            self._solver.covariance_X,
            index=self._matrix_index,
            columns=self._matrix_index,
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def covariance_Y(self) -> pd.DataFrame | None:
        """Return the covariance matrix of Y."""
        if self._solver.covariance_Y is None:
            return

        return pd.DataFrame(
            self._solver.covariance_Y,
            index=self._observation_index,
            columns=self._observation_index,
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def covariance_r(self) -> pd.DataFrame | None:
        """Return the covariance matrix of residuals."""
        if self._solver.covariance_r is None:
            return

        return pd.DataFrame(
            self._solver.covariance_r,
            index=self._observation_index,
            columns=self._observation_index,
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def increments(self) -> pd.Series | None:
        """Return increments."""
        if self._solver.increments is None:
            return

        return pd.Series(
            self._solver.increments, index=self._matrix_index, name="increments"
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def coord_increments(self) -> pd.Series | None:
        """Return fitered for just coordinate increments."""
        if self._solver.coord_increments is None:
            return

        return pd.Series(
            self._solver.coord_increments,
            index=self._matrix_coordinate_index,
            name="coord_increments",
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def increment_matrix(self) -> pd.DataFrame | None:
        """Return increment matrix."""
        if self._solver.increment_matrix is None:
            return

        return pd.DataFrame(
            self._solver.increment_matrix,
            index=self._solver.dataset.controls.index,
            columns=self._solver.dataset.controls.coordinate_columns,
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def coordinate_weights(self) -> pd.DataFrame | None:
        """Return the point weights."""
        if self._solver.coordinate_weights is None:
            return

        df = pd.Series(
            self._solver.coordinate_weights,
            index=self._matrix_coordinate_index,
            name="coordinate_weights",
        ).unstack("column", sort=False)
        df.columns = [f"w{col}" for col in df.columns]
        return df

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def obs_residuals(self) -> pd.DataFrame | None:
        """Return observation residuals."""
        if self._solver.obs_residuals is None:
            return

        df = pd.Series(
            data=self._solver.obs_residuals.flatten(),
            index=self._observation_index,
            name="obs_residuals",
        ).unstack("column", sort=False)
        df.columns = [f"v{col}" for col in df.columns]
        return df

    @property
    def n_movable_tie_points(self) -> int:
        """Return number of movable tie points."""
        return self._solver.n_movable_tie_points

    @property
    def n_fixed_tie_points(self) -> int:
        return self._solver.n_fixed_tie_points

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def residual_variance(self) -> float | None:
        return self._solver.residual_variance

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def residual_variances(self) -> pd.Series | None:
        if self._solver.residual_variances is None:
            return

        return pd.Series(
            self._solver.residual_variances,
            index=self._iter_index,
            name="residual_variance",
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def coord_correction_variance(self) -> float | None:
        return self._solver.coord_correction_variance

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def coord_correction_variances(self) -> pd.Series | None:
        if self._solver.coord_correction_variances is None:
            return

        return pd.Series(
            self._solver.coord_correction_variances,
            index=self._iter_index,
            name="coord_correction_variances",
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def residual_sigma(self) -> float | None:
        """Return value of residual sigma."""
        return self._solver.residual_sigma

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def residual_sigmas(self) -> pd.Series | None:
        if self._solver.residual_sigmas is None:
            return

        return pd.Series(
            self._solver.residual_sigmas, index=self._iter_index, name="residual_sigmas"
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def coord_correction_sigma(self) -> float | None:
        """Return value of coordinate corrections sigma."""
        return self._solver.coord_correction_sigma

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def coord_correction_sigmas(self) -> pd.Series | None:
        if self._solver.coord_correction_sigmas is None:
            return

        return pd.Series(
            self._solver.coord_correction_sigmas,
            index=pd.Index(range(1, self.n_iter + 1), name="n_iter"),
            name="coord_correction_sigmas",
        ).rename_axis("n_iter")

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def coord_corrections(self) -> pd.DataFrame:
        return pd.DataFrame(
            self._solver.coord_corrections,
            index=self.approx_coordinates.index,
            columns=[f"v{col}" for col in self.approx_coordinates.columns],
        )

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def normalized_residuals(self) -> pd.Series:
        df = pd.Series(
            self._solver.normalized_residuals,
            index=self._observation_index,
            name="normalized_obs_residuals",
        ).unstack("column", sort=False)
        df.columns = [f"v_norm {col}" for col in df.columns]
        return df

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def normalized_corrections(self) -> pd.Series:
        df = pd.Series(
            self._solver.normalized_corrections,
            index=self._matrix_coordinate_index,
        ).unstack("column", sort=False)
        df.columns = [f"v_norm {col}" for col in df.columns]
        return df

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def calculation_status(self) -> str:
        success = "Calculations succeed"
        failiture = "Calculations aborted due to SVD did not converge"
        return success if self._solver.svd_converge else failiture

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def adjusted_coordinate_sigmas(self) -> pd.DataFrame:
        return self._get_adjusted_coordinate_sigmas()

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def coordinate_error_elipses(self) -> pd.DataFrame | None:
        return self._get_coord_error_elipses()

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def adjusted_observation_values(self) -> pd.DataFrame:
        return self._get_adjusted_observation_values()

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def adjusted_observation_sigmas(self) -> pd.DataFrame:
        return self._get_adjusted_observation_sigmas()

    @refreshable_property(refresh="_is_deprecated", reset_all=True)
    def _iter_index(self):

        return pd.Index(range(1, self.n_iter + 1), name="n_iter")

    @abstractmethod
    def _get_obs_index(self) -> pd.MultiIndex:
        """Return MultiIndex object to describe observation adjustment results."""
        pass

    @abstractmethod
    def _get_matrix_coord_index(self) -> pd.MultiIndex:
        """Return MultiIndex object to describe coordinate indices in adjustment matrices."""
        pass

    @abstractmethod
    def _get_matrix_orientation_index(self) -> pd.MultiIndex | None:
        """Return MultiIndex object to describe oreientation indices in adjustment matrices."""
        pass

    @abstractmethod
    def _get_matrix_index(self) -> pd.MultiIndex:
        """Return MultiIndex object to describe adjustment matrices."""
        pass

    @abstractmethod
    def _get_adjusted_coordinate_sigmas(self) -> pd.DataFrame:
        """Return DataFrame containing adjusted coordinate sigmas."""
        pass

    @abstractmethod
    def _get_coord_error_elipses(self) -> pd.DataFrame | None:
        """Return DataFrame containing adjusted coordinate error ellipses parameters."""
        pass

    @abstractmethod
    def _get_adjusted_observation_values(self) -> pd.DataFrame:
        """Return DataFrame containing adjusted observation (measurement) values."""
        pass

    @abstractmethod
    def _get_adjusted_observation_sigmas(self) -> pd.DataFrame:
        """Return DataFrame containing adjusted observation (measurement) sigmas."""
        pass
