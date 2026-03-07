# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest

from pysurv.validators import IsNumeric
from pysurv.validators.constants import COMMENT
from pysurv.validators.utils import PySurvValidatorMode


def test_mode_from_string():
    validator = IsNumeric()
    assert validator.mode == PySurvValidatorMode.RAISE

    validator.mode = "return"
    assert validator.mode == PySurvValidatorMode.RETURN

    with pytest.raises(ValueError):
        validator.mode = "invalid"


def test_mode_from_string_on_init():
    validator = IsNumeric(mode="return")
    assert validator.mode == PySurvValidatorMode.RETURN

    with pytest.raises(ValueError):
        IsNumeric(mode="invalid")


def test_ignore_property():
    validator = IsNumeric(ignore=(1, COMMENT, lambda x: True))
    assert not validator.ignore["na"]
    assert "literals" in validator.ignore
    assert "regexes" in validator.ignore
    assert "callables" in validator.ignore
