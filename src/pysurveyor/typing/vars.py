# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Shared type-variable definitions for generic PySurveyor classes.

All bounds use forward-reference strings to avoid import cycles at runtime.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    from pysurveyor.io.abc.pysurveyor_io import PySurveyorABC


DomainClassType = TypeVar("DomainClassType", bound="PySurveyorABC")
