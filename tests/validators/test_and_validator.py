# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators import GreaterOrEqual, IsNumeric


def test_empty_data(ignore_data):
    non_negative = IsNumeric() & GreaterOrEqual(0)
    mask, values = non_negative(ignore_data.empty)

    assert mask.empty
    assert values.empty


def test_passing_output(is_numeric_data):
    with pytest.raises(ValidationError):
        non_negative = IsNumeric() & GreaterOrEqual(0)
        non_negative(is_numeric_data.convertible)

    non_negative = IsNumeric() & GreaterOrEqual(-3)
    mask, values = non_negative(is_numeric_data.convertible)

    assert mask.all()
    assert values.equals(is_numeric_data.valid)
