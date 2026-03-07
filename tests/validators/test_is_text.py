# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators import IsText


def test_is_text_empty(text_data):
    is_text = IsText()
    mask, values = is_text(text_data.empty)
    assert mask.empty
    assert mask.dtype == bool
    assert values.empty
    assert values.dtype == "str"


def test_is_text_valid(text_data):
    is_text = IsText()
    mask, values = is_text(text_data.valid)
    assert mask.all()
    assert all([validated == raw for validated, raw in zip(values, text_data.valid)])


def test_is_text_convertible(text_data):
    is_text = IsText()
    mask, values = is_text(text_data.convertible)
    assert mask.all()
    assert values.dtype == "str"


def test_is_text_with_empty(text_data):
    is_text = IsText()
    mask, values = is_text(text_data.with_empty)
    assert mask.all()
    assert values.isna().sum() == 2


def test_is_text_no_ignore(text_data):
    is_text = IsText(ignore=[])
    with pytest.raises(ValidationError):
        is_text(text_data.with_empty)


def test_is_text_ignore_literal(text_data):
    is_text = IsText(ignore=("a", "b"))
    mask, values = is_text(text_data.valid)
    assert mask.all()
    assert not values.isin(["a", "b"]).any()


def test_is_text_ignore_regex(text_data):
    is_text = IsText()
    mask, values = is_text(text_data.regex)
    assert mask.all()
    assert not values.isin(["# comment", "  # leading space", "	# leading tab"]).any()


def test_is_text_ignore_callable(text_data):
    is_text = IsText(ignore=lambda x: str(x).isdigit())
    mask, values = is_text(text_data.callable)
    assert mask.all()
    assert not values.isin(["1", "2", "3"]).any()
