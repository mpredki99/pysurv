import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators import GreaterOrEqual, IsNumeric


def test_non_negative_empty(ge_data):
    non_negative = IsNumeric() & GreaterOrEqual()
    mask, values = non_negative(ge_data.empty)

    assert mask.empty
    assert values.empty


def test_non_negative_valid(ge_data):
    non_negative = IsNumeric() & GreaterOrEqual()
    mask, values = non_negative(ge_data.valid)

    assert mask.all()
    assert values.equals(ge_data.valid)


def test_non_negative_invalid(is_numeric_data, ge_data):
    with pytest.raises(ValidationError):
        non_negative = IsNumeric() & GreaterOrEqual()
        non_negative(is_numeric_data.invalid)

    with pytest.raises(ValidationError):
        non_negative = IsNumeric() & GreaterOrEqual()
        non_negative(ge_data.invalid)


def test_non_negative_convertible(is_numeric_data):
    with pytest.raises(ValidationError):
        non_negative = IsNumeric() & GreaterOrEqual()
        non_negative(is_numeric_data.convertible)

    non_negative = IsNumeric() & GreaterOrEqual(-3)
    mask, values = non_negative(is_numeric_data.convertible)

    assert mask.all()
    assert values.equals(is_numeric_data.valid)


def test_non_negative_with_empty(is_numeric_data):
    with pytest.raises(ValidationError):
        non_negative = IsNumeric() & GreaterOrEqual()
        non_negative(is_numeric_data.with_empty)

    non_negative = IsNumeric() & GreaterOrEqual(-3)
    mask, values = non_negative(is_numeric_data.with_empty)

    assert mask.all()
    assert values.isna().equals(is_numeric_data.with_empty.isna())
