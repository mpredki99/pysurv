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
    validator = IsNumeric(ignore=(1, COMMENT, lambda x: True))
    ignored = validator.ignore_dict

    assert ignored["na"] is False
    assert ignored["literals"] == {1}
    assert ignored["regexes"] == [COMMENT]
    assert len(ignored["callables"]) == 1

    # ----------------------------------------------------------------------------------
    validator.ignore = None
    ignored = validator.ignore_dict

    assert ignored["na"] is True
    assert "literals" not in ignored
    assert "regexes" not in ignored
    assert "callables" not in ignored


# ======================================================================================
def test_contains():
    validator = Greater(0)

    assert Greater in validator
    assert Greater(1) in validator
    assert GreaterOrEqual not in validator
    assert GreaterOrEqual(0) not in validator


# ======================================================================================
def test_and():
    validator = IsNumeric() & Greater(0)

    assert isinstance(validator, AndValidator)

    # ----------------------------------------------------------------------------------
    with pytest.raises(TypeError, match="non-Validator"):
        IsNumeric() & object()


# ======================================================================================
def test_chain():
    left = IsNumeric()
    middle = Less(10)
    right = Greater(0)
    validator = left & middle & right

    assert list(validator) == [left, middle, right]
    assert len(validator) == 3
    assert validator[0] is left
    assert validator[1] is middle
    assert validator[-1] is right
    assert isinstance(validator[1:], AndValidator)
    assert list(validator[1:]) == [middle, right]
