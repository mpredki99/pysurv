# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from pysurv import Dataset
from pysurv.adjustment import Results, Solver
from pysurv.adjustment.adjustment_matrices import AdjustmentMatrices


def test_bool(adjustment_test_matrices: AdjustmentMatrices) -> None:
    """Test that beofre start calculations results are not generated."""
    solver = Solver(adjustment_test_matrices)
    results = Results(solver)

    assert solver.current_iter == 0
    assert not results


def test_iterate(
    adjustment_test_matrices: AdjustmentMatrices, adjustment_test_dataset: Dataset
) -> None:
    """Test that iterate method works properly."""
    solver = Solver(adjustment_test_matrices)
    results = Results(solver)

    solver.iterate()

    assert solver.current_iter == 1
    assert results
