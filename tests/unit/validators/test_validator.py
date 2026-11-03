# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Test for validator base functionalities."""

import pytest

from pysurveyor.validators import Greater, GreaterOrEqual, IsNumeric, Less
from pysurveyor.validators.abc.validator import AndValidator
from pysurveyor.validators.utils.constants import COMMENT
from pysurveyor.validators.utils.validator_mode import ValidatorMode


# ======================================================================================
def test_mode_on_init():
    """Validate the initial validator mode."""
    validator = IsNumeric(mode=ValidatorMode.RETURN)
    assert validator.mode == ValidatorMode.RETURN

    # ----------------------------------------------------------------------------------
    validator = IsNumeric(mode="return")
    assert validator.mode == ValidatorMode.RETURN

    # ----------------------------------------------------------------------------------
    with pytest.raises(ValueError):
        IsNumeric(mode="invalid")


# ======================================================================================
def test_mode_setter():
    """Validate changing the validator mode."""
    validator = IsNumeric()
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
def test_ignore_dict_groups_rule_types():
    """Group ignored values by their rule type."""
    validator = IsNumeric(ignore=(1, COMMENT, lambda x: True))
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
def test_contains():
    """Check validator membership in a validator chain."""
    validator = Greater(0)

    # ----------------------------------------------------------------------------------

    assert Greater in validator  # type: ignore
    assert Greater(1) in validator
    assert GreaterOrEqual not in validator  # type: ignore
    assert GreaterOrEqual(0) not in validator


# ======================================================================================
def test_and():
    """Combine validators with the ``&`` operator."""
    validator = IsNumeric() & Greater(0)

    # ----------------------------------------------------------------------------------

    assert isinstance(validator, AndValidator)

    with pytest.raises(TypeError, match="non-Validator"):
        IsNumeric() & object()  # type: ignore


# ======================================================================================
def test_chain():
    """Preserve order and slicing for chained validators."""
    left = IsNumeric()
    middle = Less(10)
    right = Greater(0)
    validator = left & middle & right

    # ----------------------------------------------------------------------------------

    assert list(validator) == [left, middle, right]
    assert len(validator) == 3
    assert validator[0] is left
    assert validator[1] is middle
    assert validator[-1] is right
    assert isinstance(validator[1:], AndValidator)
    assert list(validator[1:]) == [middle, right]  # type: ignore
