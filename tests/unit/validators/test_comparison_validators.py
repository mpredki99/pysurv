# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

import pytest

from pysurveyor.exceptions import ValidationError
from pysurveyor.validators import Equal, Greater, GreaterOrEqual, Less, LessOrEqual


# ======================================================================================
def test_threshold_as_float():
    """Convert comparison thresholds to floats and validate assignments."""
    greater = Greater(5)
    assert greater.threshold == 5.0

    # ----------------------------------------------------------------------------------

    greater.threshold = 10
    assert greater.threshold == 10.0

    # ----------------------------------------------------------------------------------

    greater.threshold = "15"  # type: ignore
    assert greater.threshold == 15.0

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError):
        greater.threshold = "non-numeric"  # type: ignore

    # ----------------------------------------------------------------------------------

    with pytest.raises(TypeError):
        greater.threshold = None  # type: ignore


# ======================================================================================
# Test validating equal
# ======================================================================================
class TestEqual:
    """Tests for ``Equal`` validator."""

    # ==================================================================================
    def test_empty(self, equal_data):
        """Accept an empty input without validation errors."""
        equal = Equal(0)
        mask, values = equal(equal_data.empty)

        # ------------------------------------------------------------------------------

        assert mask.empty
        assert mask.dtype == bool
        assert values.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, equal_data):
        """Accept values equal to the configured target."""
        equal = Equal(0)
        mask, values = equal(equal_data.valid)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(equal_data.valid)

    # ==================================================================================
    def test_invalid_type(self, equal_data):
        """Reject values with an incompatible type."""
        equal = Equal(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match="!="):
            equal(equal_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, equal_data):
        """Reject values that differ from the target."""
        equal = Equal(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match="!="):
            equal(equal_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, equal_data):
        """Allow empty values by default."""
        equal = Equal(0)
        mask, values = equal(equal_data.with_empty)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, equal_data):
        """Reject empty values when ignoring is disabled."""
        equal = Equal(0, ignore=[])

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError):
            equal(equal_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, equal_data):
        """Ignore configured literal values during equality checks."""
        equal = Equal(0, ignore=(1, 2))
        mask, values = equal(equal_data.literal)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(equal_data.literal)

    # ==================================================================================
    def test_ignore_regex(self, equal_data):
        """Ignore values matching the configured comment pattern."""
        equal = Equal(-1)
        mask, values = equal(equal_data.regex)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(equal_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, equal_data):
        """Ignore values matched by a callable."""
        equal = Equal(0, ignore=[None, lambda x: x % 2 == 1])
        mask, values = equal(equal_data.callable)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(equal_data.callable)


# ======================================================================================
# Test validating greater
# ======================================================================================
class TestGreater:
    """Tests for ``Greater`` validator."""

    # ==================================================================================
    def test_empty(self, comparison_data):
        """Accept an empty input without validation errors."""
        greater = Greater(0)
        mask, values = greater(comparison_data.empty)

        # ------------------------------------------------------------------------------

        assert mask.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, comparison_data):
        """Accept values greater than the threshold."""
        greater = Greater(-1)
        mask, values = greater(comparison_data.valid)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.valid)

    # ==================================================================================
    def test_invalid_type(self, comparison_data):
        """Reject values with an incompatible type."""
        greater = Greater(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match="Non-numeric"):
            greater(comparison_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, comparison_data):
        """Reject values at or below the threshold."""
        greater = Greater(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match="<="):
            greater(comparison_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, comparison_data):
        """Allow empty values by default."""
        greater = Greater(-1)
        mask, values = greater(comparison_data.with_empty)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, comparison_data):
        """Reject empty values when ignoring is disabled."""
        greater = Greater(-1, ignore=[])

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError):
            greater(comparison_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, comparison_data):
        """Ignore configured literal values during greater-than checks."""
        greater = Greater(0, ignore=(-1, 0))
        mask, values = greater(comparison_data.literal)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.literal)

    # ==================================================================================
    def test_ignore_regex(self, comparison_data):
        """Ignore values matching the configured comment pattern."""
        greater = Greater(-1)
        mask, values = greater(comparison_data.regex)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, comparison_data):
        """Ignore values matched by a callable."""
        greater = Greater(0, ignore=lambda x: x % 2 == 0)
        mask, values = greater(comparison_data.callable)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.callable)


# ======================================================================================
# Test validating greater or equal
# ======================================================================================
class TestGreaterOrEqual:
    """Tests for ``GreaterOrEqual`` validator."""

    # ==================================================================================
    def test_empty(self, comparison_data):
        """Accept an empty input without validation errors."""
        greater_or_equal = GreaterOrEqual(0)
        mask, values = greater_or_equal(comparison_data.empty)

        # ------------------------------------------------------------------------------

        assert mask.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, comparison_data):
        """Accept values greater than or equal to the threshold."""
        greater_or_equal = GreaterOrEqual(0)
        mask, values = greater_or_equal(comparison_data.valid)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.valid)

    # ==================================================================================
    def test_invalid_type(self, comparison_data):
        """Reject values with an incompatible type."""
        greater_or_equal = GreaterOrEqual(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match="Non-numeric"):
            greater_or_equal(comparison_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, comparison_data):
        """Reject values below the threshold."""
        greater_or_equal = GreaterOrEqual(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match="<"):
            greater_or_equal(comparison_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, comparison_data):
        """Allow empty values by default."""
        greater_or_equal = GreaterOrEqual(0)
        mask, values = greater_or_equal(comparison_data.with_empty)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, comparison_data):
        """Reject empty values when ignoring is disabled."""
        greater_or_equal = GreaterOrEqual(0, ignore=[])

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError):
            greater_or_equal(comparison_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, comparison_data):
        """Ignore configured literal values during greater-or-equal checks."""
        greater_or_equal = GreaterOrEqual(0, ignore=-1)
        mask, values = greater_or_equal(comparison_data.literal)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.literal)

    # ==================================================================================
    def test_ignore_regex(self, comparison_data):
        """Ignore values matching the configured comment pattern."""
        greater_or_equal = GreaterOrEqual(0)
        mask, values = greater_or_equal(comparison_data.regex)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, comparison_data):
        """Ignore values matched by a callable."""
        greater_or_equal = GreaterOrEqual(0, ignore=lambda x: x % 2 == 0)
        mask, values = greater_or_equal(comparison_data.callable)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.callable)


# ======================================================================================
# Test validating less
# ======================================================================================
class TestLess:
    """Tests for ``Less`` validator."""

    # ==================================================================================
    def test_empty(self, comparison_data):
        """Accept an empty input without validation errors."""
        less = Less(0)
        mask, values = less(comparison_data.empty)

        # ------------------------------------------------------------------------------

        assert mask.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, comparison_data):
        """Accept values less than the threshold."""
        less = Less(3)
        mask, values = less(comparison_data.valid)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.valid)

    # ==================================================================================
    def test_invalid_type(self, comparison_data):
        """Reject values with an incompatible type."""
        less = Less(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match="Non-numeric"):
            less(comparison_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, comparison_data):
        """Reject values at or above the threshold."""
        less = Less(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match=">="):
            less(comparison_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, comparison_data):
        """Allow empty values by default."""
        less = Less(1)
        mask, values = less(comparison_data.with_empty)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, comparison_data):
        """Reject empty values when ignoring is disabled."""
        less = Less(1, ignore=[])

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError):
            less(comparison_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, comparison_data):
        """Ignore configured literal values during less-than checks."""
        less = Less(0, ignore=(0, 1))
        mask, values = less(comparison_data.literal)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.literal)

    # ==================================================================================
    def test_ignore_regex(self, comparison_data):
        """Ignore values matching the configured comment pattern."""
        less = Less(1)
        mask, values = less(comparison_data.regex)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, comparison_data):
        """Ignore values matched by a callable."""
        less = Less(0, ignore=lambda x: x % 2 == 1)
        mask, values = less(comparison_data.callable)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.callable)


# ======================================================================================
# Test validating less or equal
# ======================================================================================
class TestLessOrEqual:
    """Tests for ``LessOrEqual`` validator."""

    # ==================================================================================
    def test_empty(self, comparison_data):
        """Accept an empty input without validation errors."""
        less_or_equal = LessOrEqual(0)
        mask, values = less_or_equal(comparison_data.empty)

        # ------------------------------------------------------------------------------

        assert mask.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, comparison_data):
        """Accept values less than or equal to the threshold."""
        less_or_equal = LessOrEqual(2)
        mask, values = less_or_equal(comparison_data.valid)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.valid)

    # ==================================================================================
    def test_invalid_type(self, comparison_data):
        """Reject values with an incompatible type."""
        less_or_equal = LessOrEqual(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match="Non-numeric"):
            less_or_equal(comparison_data.invalid_type)

    # ==================================================================================
    def test_invalid_value(self, comparison_data):
        """Reject values above the threshold."""
        less_or_equal = LessOrEqual(0)

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError, match=">"):
            less_or_equal(comparison_data.invalid_value)

    # ==================================================================================
    def test_with_empty(self, comparison_data):
        """Allow empty values by default."""
        less_or_equal = LessOrEqual(0)
        mask, values = less_or_equal(comparison_data.with_empty)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, comparison_data):
        """Reject empty values when ignoring is disabled."""
        less_or_equal = LessOrEqual(-1, ignore=[])

        # ------------------------------------------------------------------------------

        with pytest.raises(ValidationError):
            less_or_equal(comparison_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, comparison_data):
        """Ignore configured literal values during less-or-equal checks."""
        less_or_equal = LessOrEqual(0, ignore=1)
        mask, values = less_or_equal(comparison_data.literal)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.literal)

    # ==================================================================================
    def test_ignore_regex(self, comparison_data):
        """Ignore values matching the configured comment pattern."""
        less_or_equal = LessOrEqual(0)
        mask, values = less_or_equal(comparison_data.regex)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, comparison_data):
        """Ignore values matched by a callable."""
        less_or_equal = LessOrEqual(0, ignore=lambda x: x % 2 == 1)
        mask, values = less_or_equal(comparison_data.callable)

        # ------------------------------------------------------------------------------

        assert mask.all()
        assert values.equals(comparison_data.callable)
