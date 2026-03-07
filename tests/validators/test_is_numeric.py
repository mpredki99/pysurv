# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pandas as pd
import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators import IsNumeric


def test_is_numeric_empty(numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(numeric_data.empty)
    assert mask.empty
    assert mask.dtype == bool
    assert values.empty
    assert values.dtype == int


def test_is_numeric_valid(numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(numeric_data.valid)
    assert mask.all()
    assert values.equals(numeric_data.valid)


def test_is_numeric_invalid(numeric_data):
    is_numeric = IsNumeric()
    with pytest.raises(ValidationError):
        is_numeric(numeric_data.invalid)


def test_is_numeric_convertible(numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(numeric_data.convertible)
    assert mask.all()
    assert values.dtype == int


def test_is_numeric_with_empty(numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(numeric_data.with_empty)
    assert mask.all()
    assert values.isna().sum() == 4


def test_is_numeric_no_ignore(numeric_data):
    is_numeric = IsNumeric(ignore=[])
    with pytest.raises(ValidationError):
        is_numeric(numeric_data.with_empty)


def test_is_numeric_ignore_literal(numeric_data):
    validator = IsNumeric(ignore=(0, 1))
    mask, values = validator(numeric_data.valid)
    assert mask.all()
    assert not values.isin([0, 1]).any()


def test_is_numeric_ignore_regex(numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(numeric_data.regex)
    assert mask.all()
    assert not values.isin([0, 1, 2]).any()


def test_is_numeric_callable(numeric_data):
    is_numeric = IsNumeric(ignore=[pd.NA, lambda x: x % 2 == 0])
    mask, values = is_numeric(numeric_data.callable)
    assert mask.all()
    assert not values.isin([-2, 0, 2]).any()
