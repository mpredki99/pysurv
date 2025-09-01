# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from pysurv import Dataset
from pysurv.adjustment import Solver
from pysurv.adjustment.adjustment_matrices import AdjustmentMatrices


def test_iterate(
    adjustment_test_matrices: AdjustmentMatrices, adjustment_test_dataset: Dataset
) -> None:
    """Test that iterate method works properly."""
    solver = Solver(adjustment_test_matrices)

    assert solver.iterate()
    assert solver.current_iter == 1


def test_solve_observation_ordinary(
    adjustment_test_matrices: AdjustmentMatrices, adjustment_test_dataset: Dataset
) -> None:
    """Test that solve method works properly with observation ordinary method."""
    matrices = adjustment_test_matrices
    matrices.methods.obs_adj = "ordinary"
    solver = Solver(adjustment_test_matrices)

    assert solver.solve()


def test_solve_observation_weighted(
    adjustment_test_matrices: AdjustmentMatrices, adjustment_test_dataset: Dataset
) -> None:
    """Test that solve method works properly with observation weighted method."""
    matrices = adjustment_test_matrices
    matrices.methods.obs_adj = "weighted"
    solver = Solver(adjustment_test_matrices)

    assert solver.solve()


def test_solve_observation_robust(
    adjustment_test_matrices: AdjustmentMatrices, adjustment_test_dataset: Dataset
) -> None:
    """Test that solve method works properly with observation robust method."""
    matrices = adjustment_test_matrices
    matrices.methods.obs_adj = "huber"
    solver = Solver(adjustment_test_matrices)

    assert solver.solve()


def test_solve_free_adj_ordinary(
    adjustment_test_matrices: AdjustmentMatrices, adjustment_test_dataset: Dataset
) -> None:
    """Test that solve method works properly with free adjustment ordinary method."""
    matrices = adjustment_test_matrices
    matrices.methods.free_adjustment = "ordinary"
    solver = Solver(adjustment_test_matrices)

    assert solver.solve()


def test_solve_free_adj_weighted(
    adjustment_test_matrices: AdjustmentMatrices, adjustment_test_dataset: Dataset
) -> None:
    """Test that solve method works properly with free adjustment weighted method."""
    matrices = adjustment_test_matrices
    matrices.methods.free_adjustment = "weighted"
    solver = Solver(adjustment_test_matrices)

    assert solver.solve()


def test_solve_free_adj_robust(
    adjustment_test_matrices: AdjustmentMatrices, adjustment_test_dataset: Dataset
) -> None:
    """Test that solve method works properly with free adjustment robust method."""
    matrices = adjustment_test_matrices
    matrices.methods.free_adjustment = "huber"
    solver = Solver(adjustment_test_matrices)

    assert solver.solve()
