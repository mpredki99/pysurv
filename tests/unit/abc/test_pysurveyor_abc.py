# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Tests for the PySurveyor abstract-class machinery."""

from abc import abstractmethod

import pytest

from pysurveyor.abc import ABSTRACT_CLASS_ATTRIBUTE, PySurveyorABC


# ======================================================================================
# Unit Tests
# ======================================================================================
class TestPySurveyorABC:
    """Test pysurveyor.abc module custom utilities."""

    # ==================================================================================
    def test_abstract_instantiation(self):
        """Test that instantiating an abstract class raises a TypeError."""

        with pytest.raises(TypeError, match="PySurveyorABC"):
            PySurveyorABC()

    # ==================================================================================
    def test_abstract_subclass_instantiation(self):
        """Test that instantiating an abstract subclass raises a TypeError."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

        # ------------------------------------------------------------------------------

        with pytest.raises(TypeError, match="AbstractSubclass"):
            AbstractSubclass()

    # ==================================================================================
    def test_concrete_subclass_instantiation(self):
        """Test that instantiating a concrete subclass of PySurveyorABC works."""

        # ------------------------------------------------------------------------------
        # Concrete subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class ConcreteSubclass(PySurveyorABC):
            pass

        # ------------------------------------------------------------------------------

        instance = ConcreteSubclass()

        assert isinstance(instance, ConcreteSubclass)
        assert isinstance(instance, PySurveyorABC)

    # ==================================================================================
    def test_missing_abstract_attribute(self):
        """Test that a missing an abstract attribute raises a TypeError."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------

        with pytest.raises(TypeError, match=r"test_attr"):
            # --------------------------------------------------------------------------
            # Concrete subclass definition with missing abstract attribute
            # --------------------------------------------------------------------------
            class ConcreteSubclass(AbstractSubclass):
                pass

    # ==================================================================================
    def test_implemented_abstract_attribute_instantiation(self):
        """Test that a subclass with abstract attributes can be instantiated."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------
        # Concrete subclass definition with implemented abstract attribute
        # ------------------------------------------------------------------------------
        class ConcreteSubclass(AbstractSubclass):
            test_attr = "test"

        # ------------------------------------------------------------------------------

        instance = ConcreteSubclass()

        assert instance.test_attr == "test"

    # ==================================================================================
    def test_abstract_attribute_inheritance(self):
        """Test that a subclass can inherit an abstract attribute from its parent."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------
        # Abstract subclass that implements the abstract attribute
        # ------------------------------------------------------------------------------
        class ImplementedClass(AbstractSubclass):
            __pysurveyor_abstract__ = True

            test_attr = "parent-implementation"

        # ------------------------------------------------------------------------------
        # Concrete subclass with abstract attribute inherited from parent
        # ------------------------------------------------------------------------------
        class ConcreteSubclass(ImplementedClass):
            pass

        # ------------------------------------------------------------------------------

        instance = ConcreteSubclass()

        assert instance.test_attr == "parent-implementation"

    # ==================================================================================
    def test_missing_abstract_attributes_list(self):
        """Test that a TypeError lists all missing abstract attributes."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            name = ABSTRACT_CLASS_ATTRIBUTE
            version = ABSTRACT_CLASS_ATTRIBUTE
            description = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------

        with pytest.raises(TypeError, match=r"description, name, version"):
            # --------------------------------------------------------------------------
            # Concrete subclass definition with missing abstract attributes
            # --------------------------------------------------------------------------
            class ConcreteSubclass(AbstractSubclass):
                pass

    # ==================================================================================
    def test_only_missing_attribute_list(self):
        """Test that a TypeError lists only missing abstract attributes."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            name = ABSTRACT_CLASS_ATTRIBUTE
            version = ABSTRACT_CLASS_ATTRIBUTE
            description = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------

        with pytest.raises(TypeError, match=r"version"):
            # --------------------------------------------------------------------------
            # Concrete subclass definition with missing abstract attribute
            # --------------------------------------------------------------------------
            class ConcreteSubclass(AbstractSubclass):
                name = "test"
                description = "Test"

    # ==================================================================================
    @pytest.mark.parametrize("value", [None, False, 0, "", [], {}])
    def test_falsy_values_implementation(self, value: object):
        """Test that a subclass can implement an abstract attribute with falsy value."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------
        # Concrete subclass definition with missing abstract attribute
        # ------------------------------------------------------------------------------
        class ConcreteSubclass(AbstractSubclass):
            test_attr = value

        # ------------------------------------------------------------------------------

        assert ConcreteSubclass.test_attr == value

    # ==================================================================================
    def test_multiple_levels_of_inheritance(self):
        """Test that a subclass can inherit abstract attributes from multiple levels."""

        # ------------------------------------------------------------------------------
        # Abstract base subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class BaseSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            name = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------
        # Abstract intermediate subclass definition
        # ------------------------------------------------------------------------------
        class IntermediateClass(BaseSubclass):
            __pysurveyor_abstract__ = True

            version = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------
        # Concrete subclass definition with missing abstract attribute
        # ------------------------------------------------------------------------------
        class ConcreteSubclass(IntermediateClass):

            name = "test"
            version = "1.0"

        # ------------------------------------------------------------------------------

        assert ConcreteSubclass.name == "test"
        assert ConcreteSubclass.version == "1.0"

    # ==================================================================================
    def test_abstract_subclass_missing_attributes(self):
        """Test that abstract subclasses can leave abstract attributes undefined."""

        # ------------------------------------------------------------------------------
        # Abstract base subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class BaseSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------
        # Abstract intermediate subclass definition
        # ------------------------------------------------------------------------------
        class IntermediateClass(BaseSubclass):
            __pysurveyor_abstract__ = True

        # ------------------------------------------------------------------------------

        assert IntermediateClass.test_attr is ABSTRACT_CLASS_ATTRIBUTE

    # ==================================================================================
    def test_multiple_levels_of_abstract_attribute_inheritance(self):
        """Abstract attributes is enforced across multiple levels of inheritance."""

        # ------------------------------------------------------------------------------
        # Abstract base subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class BaseSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        # ------------------------------------------------------------------------------
        # Abstract intermediate subclass definition
        # ------------------------------------------------------------------------------
        class IntermediateClass(BaseSubclass):
            __pysurveyor_abstract__ = True

        # ------------------------------------------------------------------------------

        with pytest.raises(TypeError, match=r"test_attr"):
            # --------------------------------------------------------------------------
            # Concrete subclass definition with missing abstract attribute
            # --------------------------------------------------------------------------
            class ConcreteSubclass(IntermediateClass):
                pass

    # ==================================================================================
    def test_abstract_attribute_on_concrete_class(self):
        """Test that a concrete class cannot define an abstract attribute."""

        with pytest.raises(TypeError, match=r"test_attr"):
            # --------------------------------------------------------------------------
            # Concrete subclass definition with abstract attribute defined
            # --------------------------------------------------------------------------
            class ConcreteClass(PySurveyorABC):
                test_attr = ABSTRACT_CLASS_ATTRIBUTE


# ======================================================================================
# Integration Tests
# ======================================================================================
class TestIntegrationWithABCMeta:
    """Test PySurveyorABC integration with standard ABCMeta metaclass."""

    # ==================================================================================
    def test_standard_abstract_methods(self):
        """Test that a concrete subclass must implement abstract methods."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            @abstractmethod
            def survey(self) -> None: ...

        # ------------------------------------------------------------------------------
        # Concrete subclass with missing abstract method implementation
        # ------------------------------------------------------------------------------
        class ConcreteSubclass(AbstractSubclass):
            pass

        # ------------------------------------------------------------------------------

        with pytest.raises(TypeError):
            ConcreteSubclass()  # pyright: ignore[reportAbstractUsage]

    # ==================================================================================
    def test_standard_abstract_methods_implementation(self) -> None:
        """Test that a concrete subclass can implement abstract methods."""

        # ------------------------------------------------------------------------------
        # Abstract subclass of PySurveyorABC definition
        # ------------------------------------------------------------------------------
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            @abstractmethod
            def survey(self) -> None: ...

        # ------------------------------------------------------------------------------
        # Concrete subclass with abstract method implementation
        # ------------------------------------------------------------------------------
        class ConcreteSubclass(AbstractSubclass):
            def survey(self) -> None:
                pass

        # ------------------------------------------------------------------------------

        assert ConcreteSubclass() is not None
