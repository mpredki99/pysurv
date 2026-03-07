# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest

from pysurv.validators import Greater


def test_comparison_threshold_setter():
    validator = Greater(5.0)
    assert validator.threshold == 5.0

    validator.threshold = 10
    assert validator.threshold == 10.0

    validator.threshold = "15"
    assert validator.threshold == 15.0

    with pytest.raises(ValueError):
        validator.threshold = "non-numeric"

    with pytest.raises(TypeError):
        validator.threshold = None
