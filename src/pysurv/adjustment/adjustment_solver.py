# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from abc import ABC, abstractmethod

import numpy as np
import pandas as pd

from pysurv.utils.utils import reset_object_cache

from .adjustment_iteration import AdjustmentIteration
from .adjustment_matrices import AdjustmentMatrices
from .adjustment_results import AdjustmentResults
from .config_solver import config_solver


class AdjustmentSolver(ABC):
    def __init__(
        self,
        matrices: AdjustmentMatrices,
        config_solver_index: str | None = None,
        create_list_of_variances: bool = False,
    ) -> None:
        self._config_solver = self._get_config_solver(config_solver_index)
        self._matrices = matrices
        self._controls = self.dataset.controls
        self._approx_coordinates = self._controls.coordinates.copy()
        self._create_list_of_variances = create_list_of_variances

        self._n_coord_corrections = None
        self._n_movable_tie_points = self._get_n_movable_tie_points()
        self._n_fixed_tie_points = self._get_n_fixed_tie_points()
        self._residual_variances = self._get_residual_variances()
        self._coord_correction_variances = self._get_coord_correction_variances()

        self._iteration = self._get_adjustment_iteration()
        self._results = self._get_adjustment_results()

        self._matrices.methods._inject_solver(self)

    @property
    def results(self):
        """Prepare if needed and return adjustment results."""
        if self.current_iter == 0:
            return None

        if self.current_iter > self._results.n_iter:
            reset_object_cache(self._results)
        return self._results

    @property
    def matrices(self):
        """Return adjustment matrices object."""
        return self._matrices

    @property
    def methods(self):
        return self._matrices.methods

    @property
    def dataset(self):
        return self._matrices.dataset

    @property
    def approx_coordinates(self) -> pd.DataFrame:
        """Return initial approximate coordinate values."""
        return self._approx_coordinates

    @property
    def current_iter(self):
        return self._iteration.current

    @property
    def matrix_G(self):
        return self._iteration.matrix_G

    @property
    def inv_matrix_G(self):
        """Return inverse of G matrix."""
        return self._iteration.inv_matrix_G

    @property
    def cross_product(self):
        """Rerurn cross product"""
        return self._iteration.cross_product

    @property
    def covariance_X(self):
        """Return the covariance matrix of X."""
        return self._iteration.covariance_X

    @property
    def covariance_Y(self):
        """Return the covariance matrix of Y."""
        return self._iteration.covariance_Y

    @property
    def covariance_r(self):
        """Return the covariance matrix of residuals."""
        return self._iteration.covariance_r

    @property
    def increments(self):
        """Return increments."""
        return self._iteration.increments

    @property
    def coord_increments(self):
        """Return fitered for just coordinate increments."""
        return self._iteration.coord_increments

    @property
    def increment_matrix(self):
        """Return increment matrix."""
        return self._iteration.increment_matrix

    @property
    def coordinate_weights(self):
        """Return the point weights."""
        return self._iteration.coordinate_weights

    @property
    def obs_residuals(self):
        """Return observation residuals."""
        return self._iteration.obs_residuals

    @property
    def n_movable_tie_points(self) -> int:
        """Return number of movable tie points."""
        return self._n_movable_tie_points

    @property
    def n_fixed_tie_points(self) -> int:
        return self._n_fixed_tie_points

    @property
    def residual_variance(self):
        return self._iteration.residual_variance

    @property
    @abstractmethod
    def residual_variances(self):
        pass

    @property
    def n_coord_corrections(self):
        if self._n_coord_corrections is None:
            self._get_n_coord_corrections()
        return self._n_coord_corrections

    @property
    @abstractmethod
    def coord_correction_variance(self):
        pass

    @property
    @abstractmethod
    def coord_correction_variances(self):
        pass

    @property
    @abstractmethod
    def residual_sigma(self) -> np.ndarray | None:
        """Return value of residual sigma."""
        pass

    @property
    @abstractmethod
    def residual_sigmas(self) -> np.ndarray | None:
        pass

    @property
    @abstractmethod
    def coord_correction_sigma(self) -> np.ndarray | None:
        """Return value of coordinate corrections sigma."""
        pass

    @property
    @abstractmethod
    def coord_correction_sigmas(self) -> np.ndarray | None:
        pass

    @property
    @abstractmethod
    def coord_corrections(self):
        pass

    @property
    @abstractmethod
    def normalized_residuals(self) -> np.ndarray:
        pass

    @property
    @abstractmethod
    def normalized_corrections(self) -> np.ndarray:
        pass

    @property
    @abstractmethod
    def svd_converge(self):
        pass

    @abstractmethod
    def solve(self):
        """Run the adjustment process."""
        pass

    @abstractmethod
    def iterate(self):
        """Perform a single iteration of the adjustment."""
        pass

    @abstractmethod
    def update_matrices(self):
        """Update X, Y, and weight matrices."""
        pass

    @abstractmethod
    def _get_adjustment_iteration(self) -> AdjustmentIteration:
        """Returns adjustment iteration object."""
        pass

    @abstractmethod
    def _get_adjustment_results(self) -> AdjustmentResults:
        """Returns adjustment results object."""
        pass

    @abstractmethod
    def _get_residual_variances(self):
        pass

    @abstractmethod
    def _get_coord_correction_variances(self):
        pass

    @abstractmethod
    def _get_n_movable_tie_points(self) -> int:
        pass

    @abstractmethod
    def _get_n_fixed_tie_points(self) -> int:
        pass

    @abstractmethod
    def _get_n_coord_corrections(self):
        pass

    def _get_config_solver(self, index: str | None):
        """Returns config_solver row."""
        if index is None:
            index = config_solver.default_index
        return config_solver[index]
