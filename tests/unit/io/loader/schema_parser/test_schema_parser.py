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
    """Test parsing field annotations from schema headers."""

    # ==================================================================================
    def test_required_field(self, parser):
        """Mark a field as required when its header uses exclamation marks."""
        schema = parser.parse_header(["!a!"])

        assert schema._model.at["a", "required"] is True

    # ==================================================================================
    def test_optional_field(self, parser):
        """Mark a field as optional when its header uses question marks."""
        schema = parser.parse_header(["?b?"])

        assert schema._model.at["b", "required"] is False

    # ==================================================================================
    def test_unmarked_field(self, parser):
        """Leave the required setting unchanged for an unmarked field."""
        schema = parser.parse_header(["c"])

        assert schema._model.at["c", "required"] is None

    # ==================================================================================
    def test_unit(self, parser):
        """Parse and assign a unit annotation."""
        schema = parser.parse_header(["d [cm]"])

        assert schema._model.at["d", "unit"] is DistanceUnit.CENTIMETER

    # ==================================================================================
    def test_unmarked_unit(self, parser):
        """Leave the unit unset when a field has no unit annotation."""
        schema = parser.parse_header(["e"])

        assert schema._model.at["e", "unit"] is None

    # ==================================================================================
    def test_validator(self, parser):
        """Parse and assign a validator annotation."""
        schema = parser.parse_header(["a `IsNumeric()`"])

        assert isinstance(schema._model.at["a", "validator"], IsNumeric)

    # ==================================================================================
    def test_unmarked_validator(self, parser):
        """Leave the validator unset when a field has no validator annotation."""
        schema = parser.parse_header(["b"])

        assert schema._model.at["b", "validator"] is None

    # ==================================================================================
    def test_description(self, parser):
        """Parse and assign a field description."""
        schema = parser.parse_header(["c (test description)"])

        assert schema._model.at["c", "description"] == "test description"

    # ==================================================================================
    def test_unmarked_description(self, parser):
        """Leave the description unset when a field has no description annotation."""
        schema = parser.parse_header(["d"])

        assert schema._model.at["d", "description"] is None

    # ==================================================================================
    def test_all_annotations(self, parser):
        """Parse all supported annotations from a field header."""
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
    def test_all_annotations_with_spaces(self, parser):
        """Ignore surrounding whitespace when parsing all annotations."""
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
    def test_multiple_fields(self, parser):
        """Parse annotations for multiple fields in one header."""
        schema = parser.parse_header(
            ["!a!", "?b?", "c [cm]", "d `GreaterOrEqual(0)`", "e (Test description)"]
        )

        assert schema._model.at["a", "required"] is True
        assert schema._model.at["b", "required"] is False
        assert schema._model.at["c", "unit"] == DistanceUnit.CENTIMETER
        assert isinstance(schema._model.at["d", "validator"], GreaterOrEqual)
        assert schema._model.at["e", "description"] == "Test description"

    # ==================================================================================
    def test_unknown_field(self, parser):
        """Reject a header field that is absent from the schema model."""
        with pytest.raises(ValueError, match="not found"):
            parser.parse_header(["unknown"])

    # ==================================================================================
    def test_required_and_optional(self, parser):
        """Reject a field marked as both required and optional."""
        with pytest.raises(ValueError, match="both required and optional"):
            parser.parse_header(["!name! ?name?"])

    # ==================================================================================
    @pytest.mark.parametrize(
        "header",
        [
            ["!name1! !name2!"],
            ["?name1? ?name2?"],
            ["name [cm] [m]"],
            ["name `IsNumeric()` `GreaterOrEqual(0)`"],
            ["name (first description) (second description)"],
        ],
    )
    def test_duplicate_annotations(self, parser, header: list[str]):
        """Reject duplicate annotations of the same type."""
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
    def test_missing_field_name(self, parser, header: list[str]):
        """Reject headers that contain annotations but no field name."""
        with pytest.raises(ValueError):
            parser.parse_header(header)

    # ==================================================================================
    def test_no_annotations(self, parser):
        """Preserve the model when fields have no annotations."""
        schema = parser.parse_header(["a", "b", "c", "d", "e"])

        assert schema._model.equals(parser._object._model)

    # ==================================================================================
    def test_empty_header(self, parser):
        """Preserve the model when parsing an empty header."""
        schema = parser.parse_header([])

        assert schema._model.equals(parser._object._model)

    # ==================================================================================
    def test_index_named(self, parser):
        """Reject schema models whose index is not named field."""
        parser._object._model.index.name = "column"

        with pytest.raises(ValueError, match="must have index named 'field'"):
            parser.parse_header(["name"])

    # ==================================================================================
    def test_model_copy(self, parser):
        """Return a configured model copy without mutating the original."""
        result = parser.parse_header(["a", "b", "c", "d", "e"])

        pd.testing.assert_frame_equal(parser._object._model, result._model)
        assert result._model is not parser._object._model
