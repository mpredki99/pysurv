"""Tests for the SchemaParser."""

from __future__ import annotations

from copy import deepcopy

import pandas as pd
import pytest

from pysurveyor.basic.units import DistanceUnit
from pysurveyor.io.loader.schema_parser import SchemaParser
from pysurveyor.validators.comparison_validators import GreaterOrEqual
from pysurveyor.validators.conversion_validators import IsNumeric


# ======================================================================================
# Implementation of necessary classes for testing purposes
# ======================================================================================
@pytest.fixture
def parser() -> SchemaParser:
    """Create a SchemaParser with a minimal schema model."""

    class TestSchema:
        def __init__(self) -> None:
            self._model = pd.DataFrame(
                {
                    "required": [None, None, None, None, None],
                    "unit": [None, None, None, None, None],
                    "validator": [None, None, None, None, None],
                    "description": [None, None, None, None, None],
                },
                index=pd.Index(
                    ["a", "b", "c", "d", "e"],
                    name="field",
                ),
            )

        # ==============================================================================
        def copy(self, deep: bool = False) -> TestSchema:
            return deepcopy(self) if deep else self

    # ----------------------------------------------------------------------------------
    return SchemaParser(TestSchema())  # type: ignore


# ======================================================================================
# Test parse_header
# ======================================================================================
class TestParseHeader:

    # ==================================================================================
    def test_required_field(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["!a!"])

        assert schema._model.at["a", "required"] is True

    # ==================================================================================
    def test_optional_field(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["?b?"])

        assert schema._model.at["b", "required"] is False

    # ==================================================================================
    def test_unmarked_field(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["c"])

        assert schema._model.at["c", "required"] is None

    # ==================================================================================
    def test_unit(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["d [cm]"])

        assert schema._model.at["d", "unit"] is DistanceUnit.CENTIMETER

    # ==================================================================================
    def test_unmarked_unit(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["e"])

        assert schema._model.at["e", "unit"] is None

    # ==================================================================================
    def test_validator(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["a `IsNumeric()`"])

        assert isinstance(schema._model.at["a", "validator"], IsNumeric)

    # ==================================================================================
    def test_unmarked_validator(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["b"])

        assert schema._model.at["b", "validator"] is None

    # ==================================================================================
    def test_description(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["c (test description)"])

        assert schema._model.at["c", "description"] == "test description"

    # ==================================================================================
    def test_unmarked_description(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["d"])

        assert schema._model.at["d", "description"] is None

    # ==================================================================================
    def test_all_annotations(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(
            [" (Test description) !a! [m] `IsNumeric() & GreaterOrEqual(0)` "]
        )

        row = schema._model.loc["a"]

        assert row["required"] is True
        assert row["unit"] is DistanceUnit.METER
        assert isinstance(row["validator"][0], IsNumeric)
        assert isinstance(row["validator"][1], GreaterOrEqual)
        assert row["description"] == "Test description"

    # ==================================================================================
    def test_all_annotations_with_spaces(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(
            [" (  Test description  ) !a! [ m ] `  IsNumeric() & GreaterOrEqual(0)  ` "]
        )

        row = schema._model.loc["a"]

        assert row["required"] is True
        assert row["unit"] is DistanceUnit.METER
        assert isinstance(row["validator"][0], IsNumeric)
        assert isinstance(row["validator"][1], GreaterOrEqual)
        assert row["description"] == "Test description"

    # ==================================================================================
    def test_multiple_fields(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(
            ["!a!", "?b?", "c [cm]", "d `GreaterOrEqual(0)`", "e (Test description)"]
        )

        assert schema._model.at["a", "required"] is True
        assert schema._model.at["b", "required"] is False
        assert schema._model.at["c", "unit"] == DistanceUnit.CENTIMETER
        assert isinstance(schema._model.at["d", "validator"], GreaterOrEqual)
        assert schema._model.at["e", "description"] == "Test description"

    # ==================================================================================
    def test_unknown_field(self, parser: SchemaParser) -> None:
        with pytest.raises(ValueError, match="not found"):
            parser.parse_header(["unknown"])

    # ==================================================================================
    def test_required_and_optional(self, parser: SchemaParser) -> None:
        with pytest.raises(ValueError, match="both required and optional"):
            parser.parse_header(["!name! ?name?"])

    # ==================================================================================
    @pytest.mark.parametrize(
        "header",
        [
            ["!name! !age!"],
            ["?name? ?age?"],
            ["name [cm] [m]"],
            ["name `IsNumeric()` `GreaterOrEqual(0)`"],
            ["name (first description) (second description)"],
        ],
    )
    def test_duplicate_annotations(
        self, parser: SchemaParser, header: list[str]
    ) -> None:
        with pytest.raises(ValueError, match="Multiple"):
            parser.parse_header(header)

    # ==================================================================================
    @pytest.mark.parametrize(
        "header",
        [
            [""],
            ["   "],
            ["[cm]"],
            ["`IsNumeric()`"],
            ["(description)"],
        ],
    )
    def test_missing_field_name(self, parser: SchemaParser, header: list[str]) -> None:
        with pytest.raises(ValueError):
            parser.parse_header(header)

    # ==================================================================================
    def test_no_annotations(self, parser: SchemaParser) -> None:
        schema = parser.parse_header(["a", "b", "c", "d", "e"])

        assert schema._model.equals(parser._object._model)

    # ==================================================================================
    def test_empty_header(self, parser: SchemaParser) -> None:
        schema = parser.parse_header([])

        assert schema._model.equals(parser._object._model)

    # ==================================================================================
    def test_index_named(self, parser: SchemaParser) -> None:
        parser._object._model.index.name = "column"

        with pytest.raises(ValueError, match="must have index named 'field'"):
            parser.parse_header(["name"])

    # ==================================================================================
    def test_model_copy(self, parser: SchemaParser) -> None:
        result = parser.parse_header(["a", "b", "c", "d", "e"])

        pd.testing.assert_frame_equal(parser._object._model, result._model)
        assert result._model is not parser._object._model
