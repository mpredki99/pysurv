# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

import pytest

from pysurveyor.exceptions import ValidationError
from pysurveyor.validators import Between


# ======================================================================================
# Test validating between
# ======================================================================================
class TestBetween:
    """Tests for ``Between`` validator."""

    # ==================================================================================
    def test_empty(self, between_data):
        between = Between(0, 400)
        mask, values = between(between_data.empty)

        assert mask.empty
        assert mask.dtype == bool
        assert values.empty
        assert values.dtype == int

    # ==================================================================================
    def test_inclusive_both(self, between_data):
        between = Between(0, 300, inclusive="both")
        mask, values = between(between_data.valid)
        assert mask.all()
        assert values.equals(between_data.valid)

    # ==================================================================================
    def test_inclusive_left(self, between_data):
        between = Between(0, 300, inclusive="left")

        with pytest.raises(ValidationError):
            between(between_data.valid)

        # ------------------------------------------------------------------------------
        between.mode = "return"
        mask, values = between(between_data.valid)

        assert bool(mask.iloc[-1]) is False
        assert mask.iloc[:-1].all()
        assert values.equals(between_data.valid)

    # ==================================================================================
    def test_inclusive_right(self, between_data):
        between = Between(0, 300, inclusive="right")

        with pytest.raises(ValidationError):
            between(between_data.valid)

        # ------------------------------------------------------------------------------
        between.mode = "return"
        mask, values = between(between_data.valid)

        assert bool(mask.iloc[0]) is False
        assert mask.iloc[1:].all()
        assert values.equals(between_data.valid)

    # ==================================================================================
    def test_inclusive_neither(self, between_data):
        between = Between(0, 300, inclusive="neither")

        with pytest.raises(ValidationError):
            between(between_data.valid)

        # ------------------------------------------------------------------------------
        between.mode = "return"
        mask, values = between(between_data.valid)

        assert bool(mask.iloc[0]) is False
        assert bool(mask.iloc[-1]) is False
        assert mask.iloc[1:-1].all()
        assert values.equals(between_data.valid)

    # ==================================================================================
    def test_invalid_type(self, between_data):
        between = Between(0, 400)

        with pytest.raises(ValidationError, match="Non-numeric"):
            between(between_data.invalid)

    # ==================================================================================
    def test_out_of_range(self, between_data):
        between = Between(50, 150)

        with pytest.raises(ValidationError, match="out of range"):
            between(between_data.valid)

    # ==================================================================================
    def test_with_empty(self, between_data):
        between = Between(0, 400)

        mask, values = between(between_data.with_empty)
        assert mask.all()
        assert values.isna().sum() == 2

    # ==================================================================================
    def test_no_ignore(self, between_data):
        between = Between(0, 400, ignore=[])

        with pytest.raises(ValidationError):
            between(between_data.with_empty)

    # ==================================================================================
    def test_ignore_literal(self, between_data):
        between = Between(0, 200, ignore=[200, 300])

        mask, values = between(between_data.valid)
        assert mask.all()
        assert values.equals(between_data.valid)

    # ==================================================================================
    def test_ignore_regex(self, between_data):
        between = Between(0, 300)

        mask, values = between(between_data.regex)
        assert mask.all()
        assert values.equals(between_data.regex)

    # ==================================================================================
    def test_ignore_callable(self, between_data):
        between = Between(0, 110, ignore=lambda x: x > 100)

        mask, values = between(between_data.callable)
        assert mask.all()
        assert values.equals(between_data.callable)

    # ==================================================================================
    def test_inverted_bounds(self):
        with pytest.raises(ValueError, match="minimum cannot be greater"):
            Between(10, 0)

    # ==================================================================================
    def test_invalid_inclusive(self):
        with pytest.raises(ValueError, match="inclusive"):
            Between(0, 1, inclusive="middle")

    # ==================================================================================
    def test_property_setters(self, between_data):
        between = Between(0, 300, inclusive="both", mode="return")
        mask, _ = between(between_data.valid)
        assert mask.all()

        # ------------------------------------------------------------------------------
        between.minimum = 100
        between.maximum = 200
        between.inclusive = "neither"
        mask, _ = between(between_data.valid)
        assert mask.tolist() == [False, False, False, False]

        # ------------------------------------------------------------------------------
        with pytest.raises(ValueError, match="minimum cannot be greater"):
            between.minimum = 500

        with pytest.raises(ValueError, match="minimum cannot be greater"):
            between.maximum = -1

        with pytest.raises(ValueError, match="inclusive must be one of"):
            between.inclusive = "middle"
