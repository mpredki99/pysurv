# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from pysurv.adjustment.adjustment import Adjustment
from pysurv.data.dataset import Dataset

from .config import config


class Project:
    """
    Root class for managing a pysurv project, including data, configuration, and adjustment operations.

    This class provides a unified interface for handling the entire workflow from data management to
    least squares adjustment and result reporting.
    """

    def __init__(self, dataset: Dataset) -> None:
        self._config = config
        self._dataset = dataset

        self._adjustment = None

    @property
    def adjustment(self):
        """Return the adjustment object."""
        return self._adjustment

    def adjust(
        self,
        obs_adj: str = "weighted",
        obs_tuning_constants: dict | None = None,
        free_adjustment: str | None = None,
        free_adj_tuning_constants: dict | None = None,
        config_sigma_index: str | None = None,
        matrices_build_strategy: str | None = None,
        config_solver_index: str | None = None,
        create_list_of_variances: bool = False,
    ) -> None:
        """Perform least squares adjustment."""

        self._adjustment = Adjustment(
            self._dataset,
        )
        self.adjustment.methods.obs_adj = obs_adj
        self.adjustment.methods.obs_tuning_constants = obs_tuning_constants
        self.adjustment.methods.free_adjustment = free_adjustment
        self.adjustment.methods.free_adj_tuning_constants = free_adj_tuning_constants

        self.adjustment.matrices.default_sigmas_index = config_sigma_index
        self.adjustment.matrices.build_strategy = matrices_build_strategy

        self.adjustment.solver.config_index = config_solver_index
        self.adjustment.solver.create_list_of_variances = create_list_of_variances

        self._adjustment.solver.solve()
        return self._adjustment.report
