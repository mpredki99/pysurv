# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Test for validator base functionalities."""

from itertools import product

import pytest

from pysurveyor.validators import (
    Equal,
    Greater,
    GreaterOrEqual,
    IsBoolean,
    IsNumeric,
    IsText,
    Less,
    LessOrEqual,
)
from pysurveyor.validators.abc.validator import AndValidator
from pysurveyor.validators.utils.constants import COMMENT
from pysurveyor.validators.utils.validator_mode import ValidatorMode


# ======================================================================================
@pytest.mark.parametrize(
    ("validator_class"),
    [
        lambda mode: Equal(0, mode=mode),
        lambda mode: Greater(0, mode=mode),
        lambda mode: GreaterOrEqual(0, mode=mode),
        lambda mode: Less(0, mode=mode),
        lambda mode: LessOrEqual(0, mode=mode),
        IsNumeric,
        IsText,
        IsBoolean,
    ],
)
def test_mode_on_init(validator_class):
    """Validate the initial validator mode."""
    validator = validator_class(mode=ValidatorMode.RETURN)
    assert validator.mode == ValidatorMode.RETURN

    # ----------------------------------------------------------------------------------
    validator = IsNumeric(mode="return")
    assert validator.mode == ValidatorMode.RETURN

    # ----------------------------------------------------------------------------------
    with pytest.raises(ValueError):
        validator_class(mode="invalid")


# ======================================================================================
@pytest.mark.parametrize(
    ("validator_class"),
    [
        lambda: Equal(0),
        lambda: Greater(0),
        lambda: GreaterOrEqual(0),
        lambda: Less(0),
        lambda: LessOrEqual(0),
        IsNumeric,
        IsText,
        IsBoolean,
    ],
)
def test_mode_setter(validator_class):
    """Validate changing the validator mode."""
    validator = validator_class()
    assert validator.mode == ValidatorMode.RAISE

    # ----------------------------------------------------------------------------------
    validator.mode = "return"
    assert validator.mode == ValidatorMode.RETURN

    # ----------------------------------------------------------------------------------
    validator.mode = ValidatorMode.RAISE
    assert validator.mode == ValidatorMode.RAISE

    # ----------------------------------------------------------------------------------
    with pytest.raises(ValueError):
        validator.mode = "invalid"


# ======================================================================================
@pytest.mark.parametrize(
    ("validator_class"),
    [
        lambda ignore: Equal(0, ignore=ignore),
        lambda ignore: Greater(0, ignore=ignore),
        lambda ignore: GreaterOrEqual(0, ignore=ignore),
        lambda ignore: Less(0, ignore=ignore),
        lambda ignore: LessOrEqual(0, ignore=ignore),
        IsNumeric,
        IsText,
        IsBoolean,
    ],
)
def test_ignore_dict_groups_rule_types(validator_class):
    """Group ignored values by their rule type."""
    validator = validator_class(ignore=(1, COMMENT, lambda x: True))
    ignored = validator.ignore_dict

    assert ignored["empty"] is False
    assert ignored["literals"] == {1}
    assert ignored["regexes"] == [COMMENT]
    assert len(ignored["callables"]) == 1

    # ----------------------------------------------------------------------------------
    validator.ignore = None
    ignored = validator.ignore_dict

    assert ignored["empty"] is True
    assert "literals" not in ignored
    assert "regexes" not in ignored
    assert "callables" not in ignored


# ======================================================================================
@pytest.mark.parametrize(
    ("validator_class"),
    [
        lambda: Equal(0),
        lambda: Greater(0),
        lambda: GreaterOrEqual(0),
        lambda: Less(0),
        lambda: LessOrEqual(0),
        IsNumeric,
        IsText,
        IsBoolean,
    ],
)
def test_contains(validator_class):
    """Check validator membership in a validator chain."""
    validator = validator_class()

    # ----------------------------------------------------------------------------------

    assert validator_class().__class__ in validator
    assert validator_class() in validator


# ======================================================================================
@pytest.mark.parametrize(
    "validator_class",
    [
        lambda: v1() & v2()
        for v1, v2 in product(
            [
                lambda: Equal(0),
                lambda: Greater(0),
                lambda: GreaterOrEqual(0),
                lambda: Less(0),
                lambda: LessOrEqual(0),
                IsNumeric,
                IsText,
                IsBoolean,
            ],
            repeat=2,
        )
    ],
)
def test_and(validator_class):
    """Combine validators with the ``&`` operator."""
    validator = validator_class()

    # ----------------------------------------------------------------------------------

    assert isinstance(validator, AndValidator)

    with pytest.raises(TypeError, match="non-Validator"):
        validator_class() & object()  # type: ignore


# ======================================================================================
@pytest.mark.parametrize(
    ("left_validator_class", "middle_validator_class", "right_validator_class"),
    [
        (left, middle, right)
        for left, middle, right in product(
            [
                lambda: Equal(0),
                lambda: Greater(0),
                lambda: GreaterOrEqual(0),
                lambda: Less(0),
                lambda: LessOrEqual(0),
                IsNumeric,
                IsText,
                IsBoolean,
            ],
            repeat=3,
        )
    ],
)
def test_chain(left_validator_class, middle_validator_class, right_validator_class):
    """Preserve order and slicing for chained validators."""
    left = left_validator_class()
    middle = middle_validator_class()
    right = right_validator_class()
    validator = left & middle & right

    # ----------------------------------------------------------------------------------

    assert list(validator) == [left, middle, right]
    assert len(validator) == 3
    assert validator[0] is left
    assert validator[1] is middle
    assert validator[-1] is right
    assert isinstance(validator[1:], AndValidator)
    assert list(validator[1:]) == [middle, right]
