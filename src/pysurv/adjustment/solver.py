# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.


from functools import cached_property
from warnings import warn

import numpy as np

from pysurv.utils.utils import reset_object_cache
from pysurv.warnings._warnings import InvalidVarianceWarning

from .adjustment_solver import AdjustmentSolver
from .dense_iteration import DenseIteration
from .results import Results


class Solver(AdjustmentSolver):
    """Class for solving surveying adjustment task."""

    @cached_property
    def residual_variances(self):
        if not self._residual_variances:
            return self.residual_variance
        return np.array(self._residual_variances)

    @cached_property
    def coord_correction_variance(self) -> np.ndarray:
        """Return coordinate corrections variance."""
        if not self._iteration:
            return
        return self._get_coord_corrections_variance()

    @cached_property
    def coord_correction_variances(self) -> np.ndarray:
        if not self._coord_correction_variances:
            return self.coord_correction_variance
        return np.array(self._coord_correction_variances)

    @cached_property
    def residual_sigma(self) -> np.ndarray | None:
        """Return value of residual sigma."""
        if self.residual_variance is None:
            return
        return np.sqrt(self.residual_variance)

    @cached_property
    def residual_sigmas(self) -> np.ndarray | None:
        if self.residual_variances is None:
            return None
        return np.sqrt(self.residual_variances)

    @cached_property
    def coord_correction_sigma(self) -> float | None:
        """Return value of coordinate corrections sigma."""
        if self.coord_correction_variance is None:
            return
        return np.sqrt(self.coord_correction_variance)

    @cached_property
    def coord_correction_sigmas(self) -> np.ndarray | None:
        if self.coord_correction_variances is None:
            return
        return np.sqrt(self.residual_variances)

    @cached_property
    def coord_corrections(self) -> np.ndarray:
        """Return value of coordinate corrections."""
        return (self._controls.coordinates - self._approx_coordinates).values

    @cached_property
    def normalized_residuals(self) -> np.ndarray:
        """Return normalized residuals."""
        obs_residuals = self._iteration.obs_residuals.reshape(-1)
        obs_residuals_var = self._iteration.covariance_r.diagonal()
        return self._normalize_residuals(obs_residuals, obs_residuals_var)

    @cached_property
    def normalized_corrections(self) -> np.ndarray:
        """Return normalized correction values"""
        corrections = self.coord_corrections.reshape(-1)[self._iteration._coord_idx]
        corrections_var = self._iteration.covariance_X.diagonal()[
            self._iteration._coord_idx
        ]
        normalized = np.full_like(corrections, np.nan)
        normalized[self._iteration._coord_idx] = self._normalize_residuals(
            corrections, corrections_var
        )
        return normalized

    @cached_property
    def svd_converge(self) -> bool:
        """Return SVD convergence status."""
        return self._iteration.run()

    def solve(self):
        """Run the adjustment process."""
        if not self.iterate():
            return False
        return self._check_condition()

    def iterate(self):
        """Perform a single iteration of the adjustment."""
        self._prepare_iteration()

        if not self.svd_converge:
            return False
        self._process_successful_iteration()
        return True

    def update_matrices(self):
        """Update X, Y, and weight matrices."""
        self._matrices.update_xy_matrices()
        self._update_weight_matrices()

    def _update_weight_matrices(self):
        """Update weight matrices if tuning constants are present."""
        if self.methods.obs_tuning_constants:
            self._update_w_matrix()
        if self.methods.free_adj_tuning_constants:
            self._update_sw_matrix()

    def _update_w_matrix(self):
        """Update observation weight matrix."""
        obs_residuals = self._iteration.obs_residuals.reshape(-1)
        obs_residuals_var = self._iteration.covariance_r.diagonal()

        v = self._normalize_residuals(obs_residuals, obs_residuals_var)
        self._matrices.update_w_matrix(v)

    def _update_sw_matrix(self):
        """Update control point weights matrix for free adjustment."""
        increments = self._iteration.increments.reshape(-1)
        increments_var = self._iteration.covariance_X.diagonal()

        v = self._normalize_residuals(increments, increments_var)
        self._matrices.update_sw_matrix(v)

    def _normalize_residuals(self, v: np.ndarray, var_v: np.ndarray) -> np.ndarray:
        """Return normalized residuals."""
        sv = self._calculate_sigma(var_v)
        return np.divide(v, sv, out=np.full_like(v, -np.inf), where=sv > 0)

    def _calculate_sigma(self, var_v: np.ndarray) -> np.ndarray:
        """Calculate sigma based on variance."""
        invalid = var_v[var_v < 0]
        if invalid.size > 0:
            warn(
                f"{invalid.size} negative variances occured in {self.current_iter}. iteration: {invalid}.",
                InvalidVarianceWarning,
            )
            var_v = np.clip(var_v, a_min=0, a_max=None)
        return np.sqrt(var_v)

    def _prepare_iteration(self):
        """Prepare matrices for iteration."""
        if self._iteration:
            self.update_matrices()

    def _process_successful_iteration(self):
        """Process the results of a successful iteration."""
        reset_object_cache(self)
        self._update_controls()
        if self._create_list_of_variances:
            self._append_residual_variances()
        self._refresh_tuning_constants()

    def _refresh_tuning_constants(self):
        methods_to_update = ["t", "cra"]
        if self.methods.obs_adj in methods_to_update:
            self.methods._refresh_obs_tuning_constants()
        if self.methods.free_adjustment in methods_to_update:
            self.methods._refresh_free_tuning_constants()

    def _update_controls(self):
        """Update control coordinates with increments."""
        controls = self._controls
        controls.loc[:, controls.coordinate_columns] += self._iteration.increment_matrix

    def _append_residual_variances(self):
        """Append current residual and coordinate correction variances."""
        self._residual_variances.append(self.residual_variance)
        self._coord_correction_variances.append(self.coord_correction_variance)

    def _check_condition(self):
        """Check if iteration should stop or continue."""
        if self._is_increments_within_threshold() or self._is_max_iter_exceeded():
            return True
        return self.solve()

    def _is_max_iter_exceeded(self):
        """Check if current iteration number is less than max in config."""
        return self._iteration._current >= self._config_solver.max_iter

    def _is_increments_within_threshold(self):
        """Check if all increments are less than threshold in config."""
        return all(self._iteration.coord_increments <= self._config_solver.threshold)

    def _get_coord_corrections_variance(self):
        """Get value of variance of coordinate corrections."""
        if self._n_movable_tie_points > 0:
            return self._calculate_coord_coorections_variance()
        return 1

    def _calculate_coord_coorections_variance(self):
        """Calculate variance of coordinate corrections."""
        point_weights = self._iteration.coordinate_weights
        coord_corrections = self.coord_corrections.reshape(-1)

        if point_weights is not None:
            squared_corrections = (coord_corrections**2 * point_weights).sum()
        else:
            squared_corrections = (coord_corrections**2).sum()

        return np.divide(squared_corrections, self._n_movable_tie_points)

    def _get_n_movable_tie_points(self) -> int:
        """Get number of movable reference points."""
        if self._matrices.matrix_sW is None:
            return self._count_movable_tie_points_from_indexer()
        return self._count_movable_tie_points_from_sw()

    def _count_movable_tie_points_from_indexer(self) -> int:
        """
        Count how many tie points are movable based on indexer object. If sW matrix is None
        (ordinary free adjustment), than all control points are movable tie points.
        """
        return self._matrices.indexer.coordinate_indices.max().max() + 1

    def _count_movable_tie_points_from_sw(self) -> int:
        """
        Count how many tie points are movable based on control point weight matrix.
        If sW matrix is not None, than tie control points have non-zero weights.
        Control points with zero weights are not tie points.
        """
        sW = self._matrices.matrix_sW
        return sW.diagonal()[sW.diagonal() > 0].size

    def _get_n_fixed_tie_points(self) -> int:
        """Return number of fixed tie points."""
        return self._controls.coordinates.count().sum() - self._n_movable_tie_points

    def _get_residual_variances(self):
        return []

    def _get_coord_correction_variances(self):
        return []

    def _get_adjustment_iteration(self) -> DenseIteration:
        """Returns iteration object."""
        return DenseIteration(self._matrices)

    # def _get_adjustment_results(self) -> Results:
    #     return Results(self)

    def _get_n_coord_corrections(self) -> int | None:
        if self.coord_corrections is None:
            return
        return self.coord_corrections.size
