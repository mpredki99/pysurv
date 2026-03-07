# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pandas as pd
import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators.less import Less


def test_less_empty(comparison_data):
    lt = Less(0)
    mask, values = lt(comparison_data.empty)
    assert mask.empty
    assert mask.dtype == bool
    assert values.empty
    assert values.dtype == int


def test_less_valid(comparison_data):
    lt = Less(3)
    mask, values = lt(comparison_data.valid)
    assert mask.all()
    assert values.equals(comparison_data.valid)


def test_less_invalid_type(comparison_data):
    lt = Less(0)
    with pytest.raises(ValidationError):
        lt(comparison_data.invalid_type)


def test_less_invalid_value(comparison_data):
    lt = Less(0)
    with pytest.raises(ValidationError):
        lt(comparison_data.invalid_value)


def test_less_with_empty(comparison_data):
    lt = Less(1)
    mask, values = lt(comparison_data.with_empty)
    assert mask.all()
    assert values.isna().sum() == 2


def test_le_no_ignore(comparison_data):
    lt = Less(1, ignore=[])
    with pytest.raises(ValidationError):
        lt(comparison_data.with_empty)


def test_less_ignore_literal(comparison_data):
    lt = Less(0, ignore=(0, 1))
    mask, values = lt(comparison_data.invalid_value)
    assert mask.all()
    assert values.equals(comparison_data.invalid_value)


def test_le_ignore_regex(comparison_data):
    lt = Less(1)
    mask, values = lt(comparison_data.regex)
    assert mask.all()
    assert values.equals(comparison_data.regex)


def test_less_ignore_callable(comparison_data):
    lt = Less(0, ignore=lambda x: x % 2 == 1)
    mask, values = lt(comparison_data.callable)
    assert mask.all()
    assert values.equals(comparison_data.callable)
