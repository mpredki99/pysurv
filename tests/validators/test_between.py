# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pandas as pd
import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators.between import Between


def test_between_empty(between_data):
    b = Between(0, 400)
    mask, values = b(between_data.empty)
    assert mask.empty
    assert mask.dtype == bool
    assert values.empty
    assert values.dtype == int


def test_between_inclusive_both(between_data):
    b = Between(0, 300, inclusive="both")
    mask, values = b(between_data.valid)
    assert mask.all()
    assert values.equals(between_data.valid)


def test_between_inclusive_left(between_data):
    b = Between(0, 300, inclusive="left")
    with pytest.raises(ValidationError):
        b(between_data.valid)

    b = Between(0, 300, inclusive="left", mode="return")
    mask, values = b(between_data.valid)
    assert mask.iloc[-1] == False
    assert mask.iloc[:-1].all()
    assert values.equals(between_data.valid)


def test_between_inclusive_right(between_data):
    b = Between(0, 300, inclusive="right")
    with pytest.raises(ValidationError):
        b(between_data.valid)

    b = Between(0, 300, inclusive="right", mode="return")
    mask, values = b(between_data.valid)
    assert mask.iloc[0] == False
    assert mask.iloc[1:].all()
    assert values.equals(between_data.valid)


def test_between_inclusive_neither(between_data):
    b = Between(0, 300, inclusive="neither")
    with pytest.raises(ValidationError):
        b(between_data.valid)

    b = Between(0, 300, inclusive="neither", mode="return")
    mask, values = b(between_data.valid)
    assert mask.iloc[0] == False
    assert mask.iloc[-1] == False
    assert mask.iloc[1:-1].all()
    assert values.equals(between_data.valid)


def test_between_invalid_type(between_data):
    b = Between(0, 400)
    with pytest.raises(ValidationError):
        b(between_data.invalid)


def test_between_with_empty(between_data):
    b = Between(0, 400)
    mask, values = b(between_data.with_empty)
    assert mask.all()
    assert values.isna().sum() == 2


def test_between_no_ignore(between_data):
    b = Between(0, 400, ignore=[])
    with pytest.raises(ValidationError):
        b(between_data.with_empty)


def test_between_ignore_literal(between_data):
    b = Between(0, 200, ignore=[200, 300])
    mask, values = b(between_data.valid)
    assert mask.all()
    assert values.equals(between_data.valid)


def test_between_ignore_regex(between_data):
    b = Between(0, 300)
    mask, values = b(between_data.regex)
    assert mask.all()
    assert values.equals(between_data.regex)


def test_between_ignore_callable(between_data):
    b = Between(0, 110, ignore=lambda x: x > 100)
    mask, values = b(between_data.callable)
    assert mask.all()
    assert values.equals(between_data.callable)
