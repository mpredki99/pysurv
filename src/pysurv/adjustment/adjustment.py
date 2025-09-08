# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from pysurv.data.dataset import Dataset

from .adjustment_matrices import AdjustmentMatrices
from .adjustment_method_manager import AdjustmentMethodManager
from .adjustment_report import AdjustmentReport
from .adjustment_results import AdjustmentResults
from .adjustment_solver import AdjustmentSolver
from .dense_matrices import DenseMatrices
from .method_manager import MethodManager
from .report import Report
from .results import Results
from .solver import Solver


class Adjustment:
    """
    Highest level adjustment module class for running least squares adjustment
    and show results.
    """

    def __init__(
        self,
        dataset: Dataset,
        method_manager: AdjustmentMethodManager | None = None,
        matrices: AdjustmentMatrices | None = None,
        solver: AdjustmentSolver | None = None,
        results: AdjustmentResults | None = None,
        report: AdjustmentReport | None = None,
    ) -> None:
        self._dataset = dataset
        self._method_manager = method_manager or MethodManager()
        self._matrices = matrices or DenseMatrices(self._dataset, self._method_manager)
        self._solver = solver or Solver(self._matrices)
        self._results = results or Results(self._solver)
        self._report = report or Report(self._results)

    @property
    def dataset(self):
        return self._dataset

    @property
    def methods(self):
        return self._method_manager

    @property
    def matrices(self):
        return self._matrices

    @property
    def solver(self):
        return self._solver

    @property
    def results(self):
        return self._results

    @property
    def report(self):
        """Return the adjustment report."""
        return self._report
