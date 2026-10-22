# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

import pytest

from pysurveyor.exceptions import ValidationError
from pysurveyor.validators import Equal, Greater, GreaterOrEqual, Less, LessOrEqual


# ======================================================================================
def test_threshold_as_float():
    greater = Greater(5)
    assert greater.threshold == 5.0

    greater.threshold = 10
    assert greater.threshold == 10.0

    greater.threshold = "15"
    assert greater.threshold == 15.0

    with pytest.raises(ValueError):
        greater.threshold = "non-numeric"

    with pytest.raises(TypeError):
        greater.threshold = None


# ======================================================================================
# Test validating equal
# ======================================================================================
class TestEqual:
    """Tests for ``Equal`` validator."""

    # ==================================================================================
    def test_empty(self, equal_data):
        equal = Equal(0)
        mask, values = equal(equal_data.empty)

        assert mask.empty
        assert mask.dtype == bool
        assert values.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, equal_data):
        equal = Equal(0)
        mask, values = equal(equal_data.valid)

        assert mask.all()
        assert values.equals(equal_data.valid)

    # ==================================================================================
    def test_invalid_type(self, equal_data):
        equal = Equal(0)

        with pytest.raises(ValidationError, match="!="):
            equal(equal_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, equal_data):
        equal = Equal(0)

        with pytest.raises(ValidationError, match="!="):
            equal(equal_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, equal_data):
        equal = Equal(0)
        mask, values = equal(equal_data.with_empty)

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, equal_data):
        equal = Equal(0, ignore=[])
        with pytest.raises(ValidationError):
            equal(equal_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, equal_data):
        equal = Equal(0, ignore=(1, 2))
        mask, values = equal(equal_data.invalid_value)

        assert mask.all()
        assert values.equals(equal_data.invalid_value)

    # ==================================================================================
    def test_ignore_regex(self, equal_data):
        equal = Equal(-1)
        mask, values = equal(equal_data.regex)

        assert mask.all()
        assert values.equals(equal_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, equal_data):
        equal = Equal(0, ignore=[None, lambda x: x % 2 == 1])
        mask, values = equal(equal_data.callable)

        assert mask.all()
        assert values.equals(equal_data.callable)


# ======================================================================================
# Test validating greater
# ======================================================================================
class TestGreater:
    """Tests for ``Greater`` validator."""

    # ==================================================================================
    def test_empty(self, comparison_data):
        greater = Greater(0)
        mask, values = greater(comparison_data.empty)

        assert mask.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, comparison_data):
        greater = Greater(-1)
        mask, values = greater(comparison_data.valid)

        assert mask.all()
        assert values.equals(comparison_data.valid)

    # ==================================================================================
    def test_invalid_type(self, comparison_data):
        greater = Greater(0)

        with pytest.raises(ValidationError, match="Non-numeric"):
            greater(comparison_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, comparison_data):
        greater = Greater(0)

        with pytest.raises(ValidationError, match="<="):
            greater(comparison_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, comparison_data):
        greater = Greater(-1)
        mask, values = greater(comparison_data.with_empty)

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, comparison_data):
        greater = Greater(-1, ignore=[])

        with pytest.raises(ValidationError):
            greater(comparison_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, comparison_data):
        greater = Greater(0, ignore=(-1, 0))
        mask, values = greater(comparison_data.invalid_value)

        assert mask.all()
        assert values.equals(comparison_data.invalid_value)

    # ==================================================================================
    def test_ignore_regex(self, comparison_data):
        greater = Greater(-1)
        mask, values = greater(comparison_data.regex)

        assert mask.all()
        assert values.equals(comparison_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, comparison_data):
        greater = Greater(0, ignore=lambda x: x % 2 == 0)
        mask, values = greater(comparison_data.callable)

        assert mask.all()
        assert values.equals(comparison_data.callable)


# ======================================================================================
# Test validating greater or equal
# ======================================================================================
class TestGreaterOrEqual:
    """Tests for ``GreaterOrEqual`` validator."""

    # ==================================================================================
    def test_empty(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0)
        mask, values = greater_or_equal(comparison_data.empty)

        assert mask.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0)
        mask, values = greater_or_equal(comparison_data.valid)

        assert mask.all()
        assert values.equals(comparison_data.valid)

    # ==================================================================================
    def test_invalid_type(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0)

        with pytest.raises(ValidationError, match="Non-numeric"):
            greater_or_equal(comparison_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0)

        with pytest.raises(ValidationError, match="<"):
            greater_or_equal(comparison_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0)
        mask, values = greater_or_equal(comparison_data.with_empty)

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0, ignore=[])

        with pytest.raises(ValidationError):
            greater_or_equal(comparison_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0, ignore=-1)
        mask, values = greater_or_equal(comparison_data.invalid_value)

        assert mask.all()
        assert values.equals(comparison_data.invalid_value)

    # ==================================================================================
    def test_ignore_regex(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0)
        mask, values = greater_or_equal(comparison_data.regex)

        assert mask.all()
        assert values.equals(comparison_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, comparison_data):
        greater_or_equal = GreaterOrEqual(0, ignore=lambda x: x % 2 == 0)
        mask, values = greater_or_equal(comparison_data.callable)

        assert mask.all()
        assert values.equals(comparison_data.callable)


# ======================================================================================
# Test validating less
# ======================================================================================
class TestLess:
    """Tests for ``Less`` validator."""

    # ==================================================================================
    def test_empty(self, comparison_data):
        less = Less(0)
        mask, values = less(comparison_data.empty)

        assert mask.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, comparison_data):
        less = Less(3)
        mask, values = less(comparison_data.valid)

        assert mask.all()
        assert values.equals(comparison_data.valid)

    # ==================================================================================
    def test_invalid_type(self, comparison_data):
        less = Less(0)

        with pytest.raises(ValidationError, match="Non-numeric"):
            less(comparison_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, comparison_data):
        less = Less(0)

        with pytest.raises(ValidationError, match=">="):
            less(comparison_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, comparison_data):
        less = Less(1)
        mask, values = less(comparison_data.with_empty)

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, comparison_data):
        less = Less(1, ignore=[])

        with pytest.raises(ValidationError):
            less(comparison_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, comparison_data):
        less = Less(0, ignore=(0, 1))
        mask, values = less(comparison_data.invalid_value)

        assert mask.all()
        assert values.equals(comparison_data.invalid_value)

    # ==================================================================================
    def test_ignore_regex(self, comparison_data):
        less = Less(1)
        mask, values = less(comparison_data.regex)

        assert mask.all()
        assert values.equals(comparison_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, comparison_data):
        less = Less(0, ignore=lambda x: x % 2 == 1)
        mask, values = less(comparison_data.callable)

        assert mask.all()
        assert values.equals(comparison_data.callable)


# ======================================================================================
# Test validating less or equal
# ======================================================================================
class TestLessOrEqual:
    """Tests for ``LessOrEqual`` validator."""

    # ==================================================================================
    def test_empty(self, comparison_data):
        less_or_equal = LessOrEqual(0)
        mask, values = less_or_equal(comparison_data.empty)

        assert mask.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, comparison_data):
        less_or_equal = LessOrEqual(2)
        mask, values = less_or_equal(comparison_data.valid)

        assert mask.all()
        assert values.equals(comparison_data.valid)

    # ==================================================================================
    def test_invalid_type(self, comparison_data):
        less_or_equal = LessOrEqual(0)

        with pytest.raises(ValidationError, match="Non-numeric"):
            less_or_equal(comparison_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, comparison_data):
        less_or_equal = LessOrEqual(0)

        with pytest.raises(ValidationError, match=">"):
            less_or_equal(comparison_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, comparison_data):
        less_or_equal = LessOrEqual(0)
        mask, values = less_or_equal(comparison_data.with_empty)

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, comparison_data):
        less_or_equal = LessOrEqual(-1, ignore=[])

        with pytest.raises(ValidationError):
            less_or_equal(comparison_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, comparison_data):
        less_or_equal = LessOrEqual(0, ignore=1)
        mask, values = less_or_equal(comparison_data.invalid_value)

        assert mask.all()
        assert values.equals(comparison_data.invalid_value)

    # ==================================================================================
    def test_ignore_regex(self, comparison_data):
        less_or_equal = LessOrEqual(0)
        mask, values = less_or_equal(comparison_data.regex)

        assert mask.all()
        assert values.equals(comparison_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, comparison_data):
        less_or_equal = LessOrEqual(0, ignore=lambda x: x % 2 == 1)
        mask, values = less_or_equal(comparison_data.callable)

        assert mask.all()
        assert values.equals(comparison_data.callable)
