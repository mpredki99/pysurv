# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pandas as pd
import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators import GreaterOrEqual


def test_ge_empty(comparison_data):
    ge = GreaterOrEqual(0)
    mask, values = ge(comparison_data.empty)
    assert mask.empty
    assert mask.dtype == bool
    assert values.empty
    assert values.dtype == int


def test_ge_valid(comparison_data):
    ge = GreaterOrEqual(0)
    mask, values = ge(comparison_data.valid)
    assert mask.all()
    assert values.equals(comparison_data.valid)


def test_ge_invalid_type(comparison_data):
    ge = GreaterOrEqual(0)
    with pytest.raises(ValidationError):
        ge(comparison_data.invalid_type)


def test_ge_invalid_value(comparison_data):
    ge = GreaterOrEqual(0)
    with pytest.raises(ValidationError):
        ge(comparison_data.invalid_value)


def test_ge_with_empty(comparison_data):
    ge = GreaterOrEqual(0)
    mask, values = ge(comparison_data.with_empty)
    assert mask.all()
    assert values.isna().sum() == 2


def test_ge_no_ignore(comparison_data):
    ge = GreaterOrEqual(0, ignore=[])
    with pytest.raises(ValidationError):
        ge(comparison_data.with_empty)


def test_ge_ignore_literal(comparison_data):
    ge = GreaterOrEqual(0, ignore=-1)
    mask, values = ge(comparison_data.invalid_value)
    assert mask.all()
    assert values.equals(comparison_data.invalid_value)


def test_ge_ignore_regex(comparison_data):
    ge = GreaterOrEqual(0)
    mask, values = ge(comparison_data.regex)
    assert mask.all()
    assert values.equals(comparison_data.regex)


def test_ge_ignore_callable(comparison_data):
    ge = GreaterOrEqual(0, ignore=lambda x: x % 2 == 0)
    mask, values = ge(comparison_data.callable)
    assert mask.all()
    assert values.equals(comparison_data.callable)
