import pandas as pd
import pytest

from pysurv.exceptions import ValidationError
from pysurv.validators import GreaterOrEqual


def test_ge_empty(ge_data):
    ge = GreaterOrEqual()
    mask, values = ge(ge_data.empty)

    assert mask.empty
    assert values.empty


def test_ge_valid(ge_data):
    ge = GreaterOrEqual()
    mask, values = ge(ge_data.valid)

    assert mask.all()
    assert values.equals(ge_data.valid)

    with pytest.raises(ValidationError):
        ge = GreaterOrEqual(5)
        mask, values = ge(ge_data.valid)


def test_ge_invalid(ge_data):
    with pytest.raises(ValidationError):
        ge = GreaterOrEqual()
        ge(ge_data.invalid)

    ge = GreaterOrEqual(ignore={-3, -2, -1})
    mask, values = ge(ge_data.invalid)

    assert mask.all()
    assert values.equals(ge_data.invalid)


def test_ge_with_empty(ge_data):
    ge = GreaterOrEqual(ignore={pd.NA, -3})
    mask, values = ge(ge_data.with_empty)

    assert mask.all()
    assert values.equals(ge_data.with_empty)

    ge = GreaterOrEqual(ignore=[])
    with pytest.raises(ValidationError):
        ge(ge_data.with_empty)
