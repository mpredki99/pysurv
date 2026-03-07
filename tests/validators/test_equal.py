# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pandas as pd
import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators import Equal


def test_equal_empty(equal_data):
    equal = Equal(0)
    mask, values = equal(equal_data.empty)
    assert mask.empty
    assert mask.dtype == bool
    assert values.empty
    assert values.dtype == int


def test_equal_valid(equal_data):
    equal = Equal(0)
    mask, values = equal(equal_data.valid)
    assert mask.all()
    assert values.equals(equal_data.valid)


def test_equal_invalid_type(equal_data):
    equal = Equal(0)
    with pytest.raises(ValidationError):
        equal(equal_data.invalid_type)


def test_equal_invalid_value(equal_data):
    equal = Equal(0)
    with pytest.raises(ValidationError):
        equal(equal_data.invalid_value)


def test_equal_with_empty(equal_data):
    equal = Equal(0)
    mask, values = equal(equal_data.with_empty)
    assert mask.all()
    assert values.isna().sum() == 2


def test_equal_no_ignore(equal_data):
    equal = Equal(0, ignore=[])
    with pytest.raises(ValidationError):
        equal(equal_data.with_empty)


def test_equal_ignore_literal(equal_data):
    equal = Equal(2, ignore=(0, 1))
    mask, values = equal(equal_data.invalid_value)
    assert mask.all()
    assert values.equals(equal_data.invalid_value)


def test_equal_ignore_regex(equal_data):
    equal = Equal(-1)
    mask, values = equal(equal_data.regex)
    assert mask.all()
    assert values.equals(equal_data.regex)


def test_equal_callable(equal_data):
    equal = Equal(0, ignore=[pd.NA, lambda x: x % 2 == 1])
    mask, values = equal(equal_data.callable)
    assert mask.all()
    assert values.equals(equal_data.callable)
