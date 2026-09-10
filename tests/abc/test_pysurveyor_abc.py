"""Tests for the PySurveyor abstract-class machinery."""

from abc import abstractmethod

import pytest

from pysurveyor.abc import ABSTRACT_CLASS_ATTRIBUTE, PySurveyorABC


# --------------------------------------------------------------------------------------
# Unit Tests
# --------------------------------------------------------------------------------------
class TestPySurveyorABC:
    """Test pysurveyor.abc module custom utilities."""

    # ----------------------------------------------------------------------------------
    def test_abstract_base_class_cannot_be_instantiated(self) -> None:
        with pytest.raises(
            TypeError,
        ):
            PySurveyorABC()

    # ----------------------------------------------------------------------------------
    def test_abstract_subclass_cannot_be_instantiated(self) -> None:

        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

        with pytest.raises(
            TypeError,
        ):
            AbstractSubclass()

    # ----------------------------------------------------------------------------------
    def test_concrete_subclass_can_be_instantiated(self) -> None:

        class ConcreteSurveyor(PySurveyorABC):
            pass

        instance = ConcreteSurveyor()

        assert isinstance(instance, ConcreteSurveyor)
        assert isinstance(instance, PySurveyorABC)

    # ----------------------------------------------------------------------------------
    def test_missing_abstract_attribute_raises_type_error(self) -> None:

        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        with pytest.raises(TypeError, match=r"test_attr"):

            class ConcreteClass(AbstractSubclass):
                pass

    # ----------------------------------------------------------------------------------
    def test_implemented_abstract_attribute_allows_instantiation(self) -> None:

        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        class ConcreteClass(AbstractSubclass):
            test_attr = "test"

        instance = ConcreteClass()

        assert instance.test_attr == "test"

    # ----------------------------------------------------------------------------------
    def test_abstract_attribute_can_be_inherited_from_parent(self) -> None:

        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        class ImplementedClass(AbstractSubclass):
            __pysurveyor_abstract__ = True
            test_attr = "parent-implementation"

        class ConcreteClass(ImplementedClass):
            pass

        instance = ConcreteClass()

        assert instance.test_attr == "parent-implementation"

    # ----------------------------------------------------------------------------------
    def test_all_missing_abstract_attributes_are_listed(self) -> None:

        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            name = ABSTRACT_CLASS_ATTRIBUTE
            version = ABSTRACT_CLASS_ATTRIBUTE
            description = ABSTRACT_CLASS_ATTRIBUTE

        with pytest.raises(TypeError, match=r"description, name, version"):

            class ConcreteClass(AbstractSubclass):
                pass

    # ----------------------------------------------------------------------------------
    def test_only_missing_attribute_is_listed(self) -> None:

        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            name = ABSTRACT_CLASS_ATTRIBUTE
            version = ABSTRACT_CLASS_ATTRIBUTE
            description = ABSTRACT_CLASS_ATTRIBUTE

        with pytest.raises(TypeError, match=r"version"):

            class ConcreteClass(AbstractSubclass):
                name = "test"
                description = "Test"

    # ----------------------------------------------------------------------------------
    @pytest.mark.parametrize("value", [None, False, 0, "", [], {}])
    def test_falsy_values_are_valid_implementations(self, value: object) -> None:
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        class ConcreteSubclass(AbstractSubclass):
            test_attr = value

        assert ConcreteSubclass.test_attr == value

    # ----------------------------------------------------------------------------------
    def test_multiple_levels_of_inheritance_are_supported(self) -> None:

        class BaseSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            name = ABSTRACT_CLASS_ATTRIBUTE

        class IntermediateClass(BaseSubclass):
            __pysurveyor_abstract__ = True
            version = ABSTRACT_CLASS_ATTRIBUTE

        class ConcreteClass(IntermediateClass):
            __pysurveyor_abstract__ = False
            name = "test"
            version = "1.0"

        assert ConcreteClass.name == "test"
        assert ConcreteClass.version == "1.0"

    # ----------------------------------------------------------------------------------
    def test_abstract_subclass_is_allowed_to_have_missing_attributes(self) -> None:

        class BaseSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        class IntermediateClass(BaseSubclass):
            __pysurveyor_abstract__ = True

        assert IntermediateClass.test_attr is ABSTRACT_CLASS_ATTRIBUTE

    # ----------------------------------------------------------------------------------
    def test_multiple_levels_of_abstract_attribute_inheritance(self) -> None:

        class BaseSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True
            test_attr = ABSTRACT_CLASS_ATTRIBUTE

        class IntermediateClass(BaseSubclass):
            __pysurveyor_abstract__ = True

        with pytest.raises(
            TypeError,
            match=r"test_attr",
        ):

            class ConcreteSurveyor(IntermediateClass):
                __pysurveyor_abstract__ = False

    # ----------------------------------------------------------------------------------
    def test_abstract_attribute_on_concrete_class_raises_type_error(self) -> None:
        with pytest.raises(
            TypeError,
            match=r"test_attr",
        ):

            class ConcreteClass(PySurveyorABC):
                test_attr = ABSTRACT_CLASS_ATTRIBUTE


# --------------------------------------------------------------------------------------
# Integration Tests
# --------------------------------------------------------------------------------------
class TestIntegrationWithABCMeta:
    """Test PySurveyorABC integration with standard ABCMeta metaclass."""

    # ----------------------------------------------------------------------------------
    def test_standard_abstract_methods_still_prevent_instantiation(self) -> None:
        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            @abstractmethod
            def survey(self) -> None: ...

        class ConcreteClass(AbstractSubclass):
            pass

        with pytest.raises(TypeError):
            ConcreteClass()  # pyright: ignore[reportAbstractUsage]

    # ----------------------------------------------------------------------------------
    def test_standard_abstract_methods_still_work(self) -> None:

        class AbstractSubclass(PySurveyorABC):
            __pysurveyor_abstract__ = True

            @abstractmethod
            def survey(self) -> None: ...

        class ConcreteClass(AbstractSubclass):
            def survey(self) -> None:
                pass

        ConcreteClass()
