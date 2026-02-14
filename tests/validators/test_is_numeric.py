import pytest
import re
import pandas as pd

from pysurv.exceptions import ValidationError
from pysurv.validators import IsNumeric


def test_is_numeric_empty(is_numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(is_numeric_data.empty)

    assert mask.empty
    assert values.empty


def test_is_numeric_valid(is_numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(is_numeric_data.valid)

    assert mask.all()
    assert values.equals(is_numeric_data.valid)


def test_is_numeric_invalid(is_numeric_data):
    with pytest.raises(ValidationError):
        is_numeric = IsNumeric()
        is_numeric(is_numeric_data.invalid)

    is_numeric = IsNumeric(ignore=is_numeric_data.invalid)
    mask, values = is_numeric(is_numeric_data.invalid)

    assert mask.all()
    assert values.isna().all()


def test_is_numeric_convertible(is_numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(is_numeric_data.convertible)

    assert mask.all()
    assert values.equals(is_numeric_data.valid)


def test_is_numeric_with_empty(is_numeric_data):
    is_numeric = IsNumeric()
    mask, values = is_numeric(is_numeric_data.with_empty)

    assert mask.all()
    assert values.isna().equals(is_numeric_data.with_empty.isna())

    is_numeric = IsNumeric(ignore=[])
    with pytest.raises(ValidationError):
        is_numeric(is_numeric_data.with_empty)


def test_is_numeric_regex(is_numeric_data):
    is_numeric = IsNumeric(ignore=[pd.NA, re.compile(r"^[\s\t]*#")])
    mask, values = is_numeric(is_numeric_data.regex_callable)

    assert mask.all()
    assert not values.isin([0, 1]).any()


def test_is_numeric_callable(is_numeric_data):
    is_numeric = IsNumeric(ignore=[pd.NA, lambda x: isinstance(x, str)])
    mask, values = is_numeric(is_numeric_data.regex_callable)

    assert mask.all()
    assert not values.isin([0, 1, 2]).any()
