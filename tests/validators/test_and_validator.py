# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators import Greater, GreaterOrEqual, IsNumeric, Less
from pysurv.validators.pysurv_validator import AndValidator


def test_empty_data(numeric_data):
    non_negative = IsNumeric() & GreaterOrEqual(0)
    mask, values = non_negative(numeric_data.empty)

    assert mask.empty
    assert values.empty


def test_passing_output(numeric_data):
    with pytest.raises(ValidationError):
        non_negative = IsNumeric() & GreaterOrEqual(0)
        non_negative(numeric_data.convertible)

    non_negative = IsNumeric() & GreaterOrEqual(-3)
    mask, values = non_negative(numeric_data.convertible)

    assert mask.all()
    assert values.equals(numeric_data.valid)


def test_contains():
    validator = IsNumeric() & Greater(0)
    assert Greater(0) in validator
    assert GreaterOrEqual(0) not in validator

    assert Greater in validator
    assert GreaterOrEqual not in validator


def test_iter():
    left = IsNumeric()
    middle = Less(10)
    right = Greater(0)

    validator = left & middle & right

    for v, ctrl in zip(validator, [left, middle, right]):
        assert v is ctrl


def test_len():
    validator = IsNumeric() & Less(10) & Greater(0)
    assert len(validator) == 3


def test_index_int():
    validator = IsNumeric() & Less(10) & Greater(0)
    assert isinstance(validator[0], IsNumeric)
    assert isinstance(validator[1], Less)
    assert isinstance(validator[-1], Greater)


def test_index_slice():
    validator = IsNumeric() & Less(10) & Greater(0)
    assert isinstance(validator[1:], AndValidator)
