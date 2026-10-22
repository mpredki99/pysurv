# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

import pytest

from pysurveyor.exceptions import ValidationError
from pysurveyor.validators import IsNumeric, IsText


# ======================================================================================
# Test validating numeric
# ======================================================================================
class TestIsNumeric:
    """Tests for ``IsNumeric`` validator."""

    # ==================================================================================
    def test_empty(self, numeric_data):
        is_numeric = IsNumeric()
        mask, values = is_numeric(numeric_data.empty)

        assert mask.empty
        assert mask.dtype == bool
        assert values.empty
        assert values.dtype == int

    # ==================================================================================
    def test_valid(self, numeric_data):
        is_numeric = IsNumeric()
        mask, values = is_numeric(numeric_data.valid)

        assert mask.all()
        assert values.equals(numeric_data.valid)

    # ==================================================================================
    def test_invalid(self, numeric_data):
        is_numeric = IsNumeric()

        with pytest.raises(ValidationError, match="IsNumeric"):
            is_numeric(numeric_data.invalid)

    # ==================================================================================
    def test_convertible(self, numeric_data):
        is_numeric = IsNumeric()
        mask, values = is_numeric(numeric_data.convertible)

        assert mask.all()
        assert values.equals(numeric_data.valid)

    # ==================================================================================
    def test_with_empty(self, numeric_data):
        is_numeric = IsNumeric()
        mask, values = is_numeric(numeric_data.with_empty)

        assert mask.all()
        assert values.isna().sum() == 3

    # ==================================================================================
    def test_no_ignore(self, numeric_data):
        is_numeric = IsNumeric(ignore=[])

        with pytest.raises(ValidationError):
            is_numeric(numeric_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, numeric_data):
        is_numeric = IsNumeric(ignore=(0, 1))
        mask, values = is_numeric(numeric_data.valid)

        assert mask.all()
        assert not values.isin([0, 1]).any()

    # ==================================================================================
    def test_ignore_regex(self, numeric_data):
        is_numeric = IsNumeric()
        mask, values = is_numeric(numeric_data.regex)

        assert mask.all()
        assert not values.isin([0, 1, 2]).any()

    # ==================================================================================
    def test_ignore_callable(self, numeric_data):
        is_numeric = IsNumeric(ignore=[None, lambda x: x % 2 == 0])
        mask, values = is_numeric(numeric_data.callable)

        assert mask.all()
        assert not values.isin([-2, 0, 2]).any()


# ======================================================================================
# Test validating text
# ======================================================================================
class TestIsText:
    """Tests for ``IsText`` validator."""

    # ==================================================================================
    def test_empty(self, text_data):
        is_text = IsText()
        mask, values = is_text(text_data.empty)

        assert mask.empty
        assert mask.dtype == bool
        assert values.empty
        assert values.dtype == "str"

    # ==================================================================================
    def test_valid(self, text_data):
        is_text = IsText()
        mask, values = is_text(text_data.valid)

        assert mask.all()
        assert list(values) == list(text_data.valid)

    # ==================================================================================
    def test_convertible(self, text_data):
        is_text = IsText()
        mask, values = is_text(text_data.convertible)
        assert mask.all()
        assert values.dtype == "str"

    # ==================================================================================
    def test_with_empty(self, text_data):
        is_text = IsText()
        mask, values = is_text(text_data.with_empty)
        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, text_data):
        is_text = IsText(ignore=[])

        with pytest.raises(ValidationError):
            is_text(text_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, text_data):
        is_text = IsText(ignore=("a", "b"))
        mask, values = is_text(text_data.valid)

        assert mask.all()
        assert not values.isin(["a", "b"]).any()

    # ==================================================================================
    def test_ignore_regex(self, text_data):
        is_text = IsText()
        mask, values = is_text(text_data.regex)

        assert mask.all()
        assert not values.isin(
            ["# comment", "  # leading space", "\t# leading tab"]
        ).any()

    # ==================================================================================
    def test_ignore_callable(self, text_data):
        is_text = IsText(ignore=lambda x: str(x).isdigit())
        mask, values = is_text(text_data.callable)

        assert mask.all()
        assert not values.isin(["1", "2", "3"]).any()
