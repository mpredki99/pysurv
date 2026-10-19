"""Tests for I/O parser accessor."""

from typing import Any

import pytest

from pysurveyor.io.accessors import ParserAccessor


# ======================================================================================
# Implementation of necessary classes for testing purposes
# ======================================================================================
class TestSchemaParser:
    def __init__(self, schema: Any) -> None:
        self.schema = schema


# ======================================================================================
class TestSchema:
    parse = ParserAccessor(TestSchemaParser)


# ======================================================================================
class TestChildSchema(TestSchema):
    pass


# ======================================================================================
class TestAbstractSchema:
    __pysurveyor_abstract__ = True
    parse = ParserAccessor(TestSchemaParser)


# ======================================================================================
class TestSchemaArgs:
    parse = ParserAccessor(TestSchemaParser)

    def __init__(self, required: object) -> None:
        self.required = required


# ======================================================================================
# Test functions
# ======================================================================================
def test_instance_cache() -> None:
    schema = TestSchema()

    first_parser = schema.parse
    second_parser = schema.parse

    assert first_parser is second_parser
    assert schema.__dict__["_schema_parser_of_parse"] is first_parser


# ======================================================================================
def test_class_access() -> None:
    first_parser = TestSchema.parse
    second_parser = TestSchema.parse

    assert first_parser is not second_parser
    assert first_parser.schema is not second_parser.schema
    assert isinstance(first_parser.schema, TestSchema)
    assert "_schema_parser_of_parse" not in TestSchema.__dict__


# ======================================================================================
def test_subclass_bind() -> None:
    parser = TestChildSchema.parse

    assert isinstance(parser.schema, TestChildSchema)


# ======================================================================================
def test_access_on_abstract_schema() -> None:
    with pytest.raises(TypeError, match="abstract schema class"):
        TestAbstractSchema.parse


# ======================================================================================
def test_access_on_invalid_schema_constructor() -> None:
    with pytest.raises(TypeError, match="could not be instantiated") as error:
        TestSchemaArgs.parse
