# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Abstract base classes and class-attribute markers for PySurveyor.

This package provides the public abstractions used to define PySurveyor
classes with required class-level attributes.
"""

from .pysurveyor_abc import ABSTRACT_CLASS_ATTRIBUTE, PySurveyorABC

__all__ = [
    "ABSTRACT_CLASS_ATTRIBUTE",
    "PySurveyorABC",
]
