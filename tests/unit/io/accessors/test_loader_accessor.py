# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests for I/O loader accessor."""

from pysurveyor.io.accessors import LoaderAccessor


# ======================================================================================
def test_class_and_instance_loader():
    """Cache one loader for class and instance access."""

    # ----------------------------------------------------------------------------------
    # Probe class definition to simulate the loader behavior for testing
    # ----------------------------------------------------------------------------------
    class LoaderProbe:
        def __init__(self, owner):
            self.owner = owner

    # ----------------------------------------------------------------------------------
    # Domain class definition with LoaderAccessor
    # ----------------------------------------------------------------------------------
    class DomainClass:
        load = LoaderAccessor(LoaderProbe)  # type: ignore

    # ----------------------------------------------------------------------------------

    first_domain_object = DomainClass()
    second_domain_object = DomainClass()

    class_loader = DomainClass.load  # type: ignore
    first_instance_loader = first_domain_object.load  # type: ignore
    second_instance_loader = second_domain_object.load  # type: ignore

    # ----------------------------------------------------------------------------------

    assert class_loader is first_instance_loader is second_instance_loader
    assert class_loader.owner is DomainClass
    assert DomainClass.__dict__["_loader_of_DomainClass"] is class_loader


# ======================================================================================
def test_subclasses_loader():
    """Bind and cache a separate loader for each subclass."""

    # ----------------------------------------------------------------------------------
    # Probe class definition to simulate the loader behavior for testing
    # ----------------------------------------------------------------------------------
    class LoaderProbe:
        def __init__(self, owner):
            self.owner = owner

    # ----------------------------------------------------------------------------------
    # Domain class definition with LoaderAccessor
    # ----------------------------------------------------------------------------------
    class DomainClass:
        load = LoaderAccessor(LoaderProbe)  # type: ignore

    # ----------------------------------------------------------------------------------
    # Child domain class definition to test loader behavior in subclasses
    # ----------------------------------------------------------------------------------
    class ChildDomainClass(DomainClass):
        pass

    # ----------------------------------------------------------------------------------

    parent_loader = DomainClass.load  # type: ignore
    child_loader = ChildDomainClass.load  # type: ignore

    # ----------------------------------------------------------------------------------

    assert child_loader is not parent_loader
    assert child_loader.owner is ChildDomainClass
    assert ChildDomainClass.__dict__["_loader_of_ChildDomainClass"] is child_loader
