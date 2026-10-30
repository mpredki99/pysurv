# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests for pysurveyor.validators.parser."""

import re

import pytest

from pysurveyor.backend.backend import NA
from pysurveyor.exceptions import ParseError
from pysurveyor.validators import (
    COMMENT,
    Between,
    Equal,
    Greater,
    GreaterOrEqual,
    IsBoolean,
    IsNumeric,
    IsText,
    Less,
    LessOrEqual,
    parse_validator,
)
from pysurveyor.validators.abc.validator import Validator


# ======================================================================================
# Basic construction and return type
# ======================================================================================
@pytest.mark.parametrize(
    ("expression", "expected_type"),
    [
        ("IsNumeric()", IsNumeric),
        ("IsBoolean()", IsBoolean),
        ("IsText()", IsText),
        ("Equal(3)", Equal),
        ("Greater(3)", Greater),
        ("GreaterOrEqual(3)", GreaterOrEqual),
        ("Less(3)", Less),
        ("LessOrEqual(3)", LessOrEqual),
        ("Between(1, 15)", Between),
    ],
)
def test_constructs_validators(expression, expected_type):
    result = parse_validator(expression)
    assert isinstance(result, expected_type)


# ======================================================================================
def test_between_arguments_are_parsed():
    result = parse_validator("Between(1, 15)")
    assert result.minimum == 1.0  # type: ignore
    assert result.maximum == 15.0  # type: ignore


# ======================================================================================
def test_comparison_threshold_is_parsed():
    result = parse_validator("GreaterOrEqual(5)")
    assert result.threshold == 5.0  # type: ignore


# ======================================================================================
# ignore: scalar, iterable, strings, callables and regex
# ======================================================================================
@pytest.mark.parametrize(
    "expression",
    [
        "IsNumeric(ignore=None)",
        "IsNumeric(ignore=[])",
        "IsNumeric(ignore=5)",
        "IsNumeric(ignore='N/A')",
        "IsNumeric(ignore=b'NA')",
        "IsNumeric(ignore=[None, 5, COMMENT])",
        "IsNumeric(ignore=(None, 5, 'N/A'))",
        "IsNumeric(ignore={1, 2, 3})",
    ],
)
def test_supported_ignore_scalar_and_iterable(expression):
    result = parse_validator(expression)
    assert isinstance(result, IsNumeric)


# ======================================================================================
def test_ignore_string():
    result = parse_validator("IsNumeric(ignore='N/A')")
    assert result.ignore == "N/A"


# ======================================================================================
def test_ignore_comment():
    result = parse_validator("IsNumeric(ignore=COMMENT)")
    assert result.ignore is COMMENT


# ======================================================================================
def test_ignore_lambda_callable():
    result = parse_validator("IsNumeric(ignore=lambda x: x % 2 == 0)")

    assert callable(result.ignore)
    assert result.ignore(4) is True
    assert result.ignore(5) is False


# ======================================================================================
def test_ignore_regex():
    result = parse_validator("IsNumeric(ignore=re.compile(r'^N/A$', re.IGNORECASE))")
    assert isinstance(result.ignore, re.Pattern)
    assert result.ignore.match("n/a")


# ======================================================================================
@pytest.mark.parametrize(
    ("expression", "expected_ignore_empty"),
    [
        ("IsNumeric(ignore=[])", False),
        ("IsNumeric(ignore=pd.NA)", True),
        ("IsNumeric(ignore=pandas.NA)", True),
        ("IsNumeric(ignore=None)", True),
        ("IsNumeric(ignore=[1, pd.NA])", True),
        ("IsNumeric(ignore=[1, pandas.NA])", True),
        ("IsNumeric(ignore=[1, None])", True),
    ],
)
def test_ignore_empty(expression, expected_ignore_empty):
    result = parse_validator(expression)
    assert result.ignore_dict["empty"] is expected_ignore_empty


# ======================================================================================
def test_lambda_unknown_variable():
    with pytest.raises(ParseError, match="Unknown name"):
        parse_validator("IsNumeric(ignore=lambda x: y > 1)")


# ======================================================================================
def test_and_validator():
    result = parse_validator("IsNumeric() & GreaterOrEqual(5) & LessOrEqual(15)")
    assert isinstance(result, Validator)
    assert len(result) == 3
    assert isinstance(result[0], IsNumeric)
    assert isinstance(result[1], GreaterOrEqual)
    assert isinstance(result[2], LessOrEqual)


# ======================================================================================
@pytest.mark.parametrize(
    "expression",
    [
        "",
        "   ",
        "IsNumeric(",
        "NotAValidator()",
        "IsNumeric() + Greater(1)",
        "IsNumeric() | Greater(1)",
        "IsNumeric(ignore=)",
        "pd.NA",
        "re.compile('x')",
        "IsNumeric(ignore={**{'x': 1} })",
        "Between(1)",
        "GreaterOrEqual(1, 2)",
        "IsNumeric(unknown_option=True)",
    ],
)
def test_invalid_expressions(expression):
    with pytest.raises(ParseError):
        parse_validator(expression)
