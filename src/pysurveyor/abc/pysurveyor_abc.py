from __future__ import annotations

from abc import ABC
from typing import ClassVar, Self


# --------------------------------------------------------------------------------------
# Abstract base class for PySurveyor framework
# --------------------------------------------------------------------------------------
class PySurveyorABC(ABC):
    """Base class for PySurveyor classes that declare abstract class attributes.

    ``PySurveyorABC`` extends :class:`abc.ABC` with support for abstract
    *class attributes*. Python's standard ``abc`` machinery handles abstract
    methods and properties, but does not provide an equivalent mechanism for
    requiring subclasses to define class-level configuration attributes.

    A subclass is considered abstract when its ``__pysurveyor_abstract__``
    class attribute is ``True``. If sublclass do not define this attribute at all,
    it is considered to be concrete. Abstract subclasses are allowed to inherit
    or declare additional abstract class attributes without implementing them.

    A concrete subclass must provide a value for every class attribute that
    is declared as :data:`ABSTRACT_CLASS_ATTRIBUTE` somewhere in its method
    resolution order (MRO). If any such attribute remains unresolved,
    :class:`TypeError` is raised when the subclass is created.

    The base class itself, and any subclass explicitly marked as abstract,
    cannot be instantiated.
    """

    __pysurveyor_abstract__: ClassVar[bool] = True

    # ----------------------------------------------------------------------------------
    def __new__(cls, *args: object, **kwargs: object) -> Self:
        """Create an instance of a concrete PySurveyor class.

        Parameters
        ----------
        *args
            Positional arguments passed to the class constructor.
        **kwargs
            Keyword arguments passed to the class constructor.

        Returns
        -------
        PySurveyorABC
            A newly allocated instance of ``cls``.

        Raises
        ------
        TypeError
            If ``cls`` is marked as abstract.

        Notes
        -----
        Abstract class attributes are validated when a concrete subclass is
        created, rather than when an instance is created. This check remains
        here to prevent direct instantiation of classes explicitly marked as
        abstract.
        """
        if cls.__dict__.get("__pysurveyor_abstract__", False):
            raise TypeError(f"Can't instantiate abstract class {cls.__name__}")

        return super().__new__(cls)

    # ----------------------------------------------------------------------------------
    def __init_subclass__(cls, **kwargs: object) -> None:
        """Validate abstract class attributes when a subclass is created.

        Parameters
        ----------
        **kwargs
            Keyword arguments forwarded to :meth:`object.__init_subclass__`.

        Raises
        ------
        TypeError
            If ``cls`` is concrete but does not implement abstract class
            attributes inherited from its base classes.

        Notes
        -----
        Validation is performed at class-definition time. This provides an
        immediate and deterministic failure when a concrete subclass is
        incomplete, instead of deferring the error until the class is
        instantiated.
        """
        super().__init_subclass__(**kwargs)

        if cls.__dict__.get("__pysurveyor_abstract__", False):
            return

        abstract_attributes = _get_abstract_class_attributes(cls)

        missing_attributes = {
            attribute_name
            for attribute_name in abstract_attributes
            if getattr(cls, attribute_name, ABSTRACT_CLASS_ATTRIBUTE)
            is ABSTRACT_CLASS_ATTRIBUTE
        }

        if missing_attributes:
            attributes = ", ".join(sorted(missing_attributes))
            raise TypeError(
                f"Can't create concrete class {cls.__name__} with "
                f"unresolved abstract class attribute(s): {attributes}"
            )


# --------------------------------------------------------------------------------------
# Abstract class attribute utilities
# --------------------------------------------------------------------------------------
class _AbstractClassAttribute:
    """Sentinel used to mark a class attribute as abstract.

    The sentinel is intentionally a unique object rather than a special value
    such as ``None`` so that ``None`` remains a valid implementation value for
    an abstract class attribute.
    """

    def __repr__(self) -> str:
        return "<abstract class attribute>"


ABSTRACT_CLASS_ATTRIBUTE = _AbstractClassAttribute()


def _get_abstract_class_attributes(cls: type) -> set[str]:
    """Return abstract class attributes declared anywhere in ``cls``'s MRO.

    Parameters
    ----------
    cls
        Class whose method resolution order should be inspected.

    Returns
    -------
    set[str]
        Names of all class attributes whose value is
        :data:`ABSTRACT_CLASS_ATTRIBUTE`.

    Notes
    -----
    The complete MRO is inspected rather than only ``cls.__dict__`` so that
    abstract class attributes declared by any ancestor are inherited as
    requirements by concrete subclasses.
    """
    return {
        attribute_name
        for base_class in cls.__mro__
        for attribute_name, value in base_class.__dict__.items()
        if value is ABSTRACT_CLASS_ATTRIBUTE
    }
