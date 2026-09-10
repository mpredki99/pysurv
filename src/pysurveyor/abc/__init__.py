"""Abstract base classes and class-attribute markers for PySurveyor.

This package provides the public abstractions used to define PySurveyor
classes with required class-level attributes.
"""

from .pysurveyor_abc import ABSTRACT_CLASS_ATTRIBUTE, PySurveyorABC

__all__ = [
    "ABSTRACT_CLASS_ATTRIBUTE",
    "PySurveyorABC",
]
