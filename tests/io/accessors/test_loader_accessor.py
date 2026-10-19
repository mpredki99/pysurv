"""Tests for I/O loader accessor."""

from typing import Any

from pysurveyor.io.accessors import LoaderAccessor


# ======================================================================================
# Implementation of necessary classes for testing purposes
# ======================================================================================
class TestLoader:
    def __init__(self, owner: type[Any]) -> None:
        self.owner = owner


# ======================================================================================
class TestDomainClass:
    load = LoaderAccessor(TestLoader)


# ======================================================================================
class TestChildDomainClass(TestDomainClass):
    pass


# ======================================================================================
# Test functions
# ======================================================================================
def test_class_and_instance_loader() -> None:
    domain_object = TestDomainClass()

    class_loader = TestDomainClass.load
    instance_loader = domain_object.load

    assert class_loader is instance_loader
    assert class_loader.owner is TestDomainClass
    assert TestDomainClass.__dict__["_loader_of_load"] is class_loader


# ======================================================================================
def test_subclasses_loader() -> None:
    parent_loader = TestDomainClass.load
    child_loader = TestChildDomainClass.load

    assert child_loader is TestChildDomainClass.load
    assert child_loader is not parent_loader
    assert child_loader.owner is TestChildDomainClass
    assert TestChildDomainClass.__dict__["_loader_of_load"] is child_loader
