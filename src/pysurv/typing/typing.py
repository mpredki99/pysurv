# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from enum import StrEnum
from typing import Literal

AngleUnit = Literal["grad", "gon", "deg", "rad"]


class PySurvUnit(StrEnum):
    pass


class DistanceUnit(PySurvUnit):
    METERS = "meters"


class AngleUnit(PySurvUnit):
    GRAD = "grad"
    DEG = "deg"
    RAD = "rad"
