# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest

from pysurv.schema.field_compilers import compile_unit, compile_validator
from pysurv.typing import AngleUnit, DistanceUnit, PySurvUnit
from pysurv.validators import Between, IsNumeric, PySurvValidator


# --------------------------------------------------------------------------------------
# Compile PySurv Validators
# --------------------------------------------------------------------------------------
def test_compile_validator_valid():
    expr = "Between(0, 400)"
    validator = compile_validator(expr)
    assert isinstance(validator, Between)


def test_compile_validator_valid():
    expr = "NotExistingValidator(1)"
    validator = compile_validator(expr)
    assert isinstance(validator, str)
    assert validator == expr


def test_compile_validator_from_validator():
    expr = IsNumeric()
    validator = compile_validator(expr)
    assert expr is validator


def test_compile_validator_deepcopy():
    expr = "Between(0, 400)"
    validator_1 = compile_validator(expr)
    validator_2 = compile_validator(expr)
    assert isinstance(validator_1, PySurvValidator)
    assert isinstance(validator_2, PySurvValidator)
    assert validator_1 is not validator_2


# --------------------------------------------------------------------------------------
# Compile PySurv Units
# --------------------------------------------------------------------------------------
def test_compile_unit_distance():
    unit = compile_unit("meters")
    assert unit is DistanceUnit.METERS


def test_compile_unit_angle():
    unit = compile_unit("grad")
    assert unit is AngleUnit.GRAD


def test_compile_unit_from_unit():
    expr = DistanceUnit.METERS
    unit = compile_unit(expr)
    assert isinstance(unit, PySurvUnit)
    assert unit is expr


@pytest.mark.parametrize("value", [None, "", "   ", 123])
def test_compile_unit_invalid(value):
    assert compile_unit(value) == value
