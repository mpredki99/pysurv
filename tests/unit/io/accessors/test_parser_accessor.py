# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests for I/O parser accessor."""

import pytest

from pysurveyor.abc import PySurveyorABC
from pysurveyor.io.accessors import ParserAccessor


# ======================================================================================
def test_instance_cache():
    """Cache a parser on each schema instance."""

    # ----------------------------------------------------------------------------------
    # Probe class definition to simulate the parser behavior for testing
    # ----------------------------------------------------------------------------------
    class SchemaParserProbe:
        def __init__(self, schema):
            self.schema = schema

    # ----------------------------------------------------------------------------------
    # Domain class definition with ParserAccessor
    # ----------------------------------------------------------------------------------
    class DomainSchema:
        parse = ParserAccessor(SchemaParserProbe)  # type: ignore

    # ----------------------------------------------------------------------------------

    domain_schema = DomainSchema()

    first_parser = domain_schema.parse  # type: ignore
    second_parser = domain_schema.parse  # type: ignore

    # ----------------------------------------------------------------------------------

    assert first_parser is second_parser
    assert domain_schema.__dict__["_schema_parser_of_DomainSchema"] is first_parser


# ======================================================================================
def test_class_access():
    """Create a fresh parser for each class-level access."""

    # ----------------------------------------------------------------------------------
    # Probe class definition to simulate the parser behavior for testing
    # ----------------------------------------------------------------------------------
    class SchemaParserProbe:
        def __init__(self, schema):
            self.schema = schema

    # ----------------------------------------------------------------------------------
    # Domain class definition with ParserAccessor
    # ----------------------------------------------------------------------------------
    class DomainSchema:
        parse = ParserAccessor(SchemaParserProbe)  # type: ignore

    # ----------------------------------------------------------------------------------

    first_parser = DomainSchema.parse  # type: ignore
    second_parser = DomainSchema.parse  # type: ignore

    # ----------------------------------------------------------------------------------

    assert first_parser is not second_parser
    assert first_parser.schema is not second_parser.schema
    assert isinstance(first_parser.schema, DomainSchema)
    assert "_schema_parser_of_DomainSchema" not in DomainSchema.__dict__


# ======================================================================================
def test_subclass_bind():
    """Bind class-level parser access to the requested subclass."""

    # ----------------------------------------------------------------------------------
    # Probe class definition to simulate the parser behavior for testing
    # ----------------------------------------------------------------------------------
    class SchemaParserProbe:
        def __init__(self, schema):
            self.schema = schema

    # ----------------------------------------------------------------------------------
    # Domain class definition with ParserAccessor
    # ----------------------------------------------------------------------------------
    class DomainSchema:
        parse = ParserAccessor(SchemaParserProbe)  # type: ignore

    # ----------------------------------------------------------------------------------
    # Child domain class definition to test parser behavior in subclasses
    # ----------------------------------------------------------------------------------
    class ChildSchema(DomainSchema):
        pass

    # ----------------------------------------------------------------------------------

    parser = ChildSchema.parse  # type: ignore

    assert isinstance(parser.schema, ChildSchema)


# ======================================================================================
def test_access_on_abstract_schema():
    """Reject parser access on abstract schema classes."""

    # ----------------------------------------------------------------------------------
    # Probe class definition to simulate the parser behavior for testing
    # ----------------------------------------------------------------------------------
    class SchemaParserProbe:
        def __init__(self, schema):
            self.schema = schema

    # ----------------------------------------------------------------------------------
    # Abstract schema class definition with parser accessor
    # ----------------------------------------------------------------------------------
    class AbstractSchema(PySurveyorABC):
        __pysurveyor_abstract__ = True

        parse = ParserAccessor(SchemaParserProbe)  # type: ignore

    # ----------------------------------------------------------------------------------

    with pytest.raises(TypeError, match="abstract schema class"):
        AbstractSchema.parse  # type: ignore


# ======================================================================================
def test_access_on_invalid_schema_constructor() -> None:
    """Explain when class-level parser binding cannot instantiate a schema."""

    # ----------------------------------------------------------------------------------
    # Probe class definition to simulate the parser behavior for testing
    # ----------------------------------------------------------------------------------
    class SchemaParserProbe:
        def __init__(self, schema):
            self.schema = schema

    # ----------------------------------------------------------------------------------
    # Domain schema class definition with required constructor argument
    # ----------------------------------------------------------------------------------
    class DomainSchema:
        parse = ParserAccessor(SchemaParserProbe)  # type: ignore

        def __init__(self, required: object) -> None:
            self.required = required

    # ----------------------------------------------------------------------------------

    with pytest.raises(TypeError, match="could not be instantiated"):
        DomainSchema.parse  # type: ignore
