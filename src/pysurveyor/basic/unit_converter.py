# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Unit conversion utilities for angles and lengths.

This module provides an API for converting between angular and linear
units using NumPy-compatible numerical operations.

Supported angle units
---------------------
- `radian`: radians
- `milliradian`: milliradians
- `gon` and `grad`: gradians (400 gon per full revolution)
- `c`: centesimal angular minutes (1/100 gon)
- `cc`: centesimal angular seconds (1/10000 gon)
- `degree`: decimal degrees
- `arcminute`: angular minutes (1/60 degree)
- `arcsecond`: angular seconds (1/3600 degree)
- `dms`: degrees, minutes, and seconds represented by arrays of shape
  (3,) or (N, 3) - individual value is `[degrees, minutes, seconds]`
- ``mm_over_meter``: slope expressed as millimeters per meter

Supported length units
----------------------
- Metric: `kilometer`, `hectometer`, `meter`, `decimeter`,
  `centimeter`, `millimeters`, and `micrometer`
- Imperial: `mile`, `yard`, `feet`, `inch`, and `thou`

Usage
-----
Convert a scalar angle:

    convert.degrees.to.rad(180)  # 3.141592653589793

Convert a NumPy-compatible array of lengths:

    convert.kilometer.to.meter([1, 2, 3])  # array([1000., 2000., 3000.])

Convert degrees, minutes, and seconds:

    convert.dms.to.degrees([12, 30, 0])  # 12.5

Convert a slope to an angle::

    convert.mm_over_meter.to.degrees(100)  # approximately 5.7106

Return values
-------------
Scalar numeric inputs are returned as Python ``float`` values. Array-like inputs
return NumPy arrays with ``float64`` precision. DMS conversions use arrays.

Raises
------
AttributeError
    If a source or target unit is invalid, or if a conversion between
    angular and length units is attempted.

ValueError
    If an unsupported dimension is specified internally, or if DMS input
    has an invalid shape or negative minutes or seconds.
"""

from __future__ import annotations

from numbers import Number
from typing import overload

import numpy as np
from numpy.typing import ArrayLike

from pysurveyor.basic.units import AngleUnit, DistanceUnit
from pysurveyor.basic.utils import FloatArray, as_array, normalize_numeric_type

# ======================================================================================
# Unit conversion factors
# ======================================================================================
ANGLE_FACTORS: dict[AngleUnit, float] = {
    AngleUnit.RADIAN: 1.0,
    AngleUnit.MILLIRADIAN: 0.001,
    AngleUnit.GON: np.pi / 200.0,
    AngleUnit.C: np.pi / 20000.0,
    AngleUnit.CC: np.pi / 2000000.0,
    AngleUnit.DEGREE: np.pi / 180.0,
    AngleUnit.ARCMINUTE: np.pi / 10800.0,
    AngleUnit.ARCSECOND: np.pi / 648000.0,
}

# ======================================================================================
DISTANCE_FACTORS: dict[DistanceUnit, float] = {
    DistanceUnit.KILOMETER: 1000.0,
    DistanceUnit.HECTOMETER: 100.0,
    DistanceUnit.METER: 1.0,
    DistanceUnit.DECIMETER: 0.1,
    DistanceUnit.CENTIMETER: 0.01,
    DistanceUnit.MILLIMETER: 0.001,
    DistanceUnit.MICROMETER: 0.000001,
    # ------------------------------
    DistanceUnit.MILE: 1609.344,
    DistanceUnit.YARD: 0.9144,
    DistanceUnit.FEET: 0.3048,
    DistanceUnit.INCH: 0.0254,
    DistanceUnit.THOU: 0.0000254,
}

# ======================================================================================
SPECIAL_ANGLE_UNITS = {
    AngleUnit.DMS,
    AngleUnit.MM_OVER_METER,
}


# ======================================================================================
# Unit normalization
# ======================================================================================
def _resolve_unit(unit: str | AngleUnit | DistanceUnit) -> AngleUnit | DistanceUnit:
    """Resolve a unit string or enum member to its canonical enum member."""
    if isinstance(unit, (AngleUnit, DistanceUnit)):
        return unit

    if not isinstance(unit, str):
        raise TypeError(
            f"Unit must be a string or unit enum, got {type(unit).__name__}."
        )

    name = unit.lower()

    for enum_type in (AngleUnit, DistanceUnit):
        try:
            return enum_type(name)
        except ValueError:
            continue

    raise ValueError(f"Unknown unit: {unit!r}.")


# ======================================================================================
# Conversion helper functions
# ======================================================================================
def _angle_to_radians(value: object, unit: AngleUnit) -> float | FloatArray:
    """Convert an angle from the specified unit to radians."""

    if unit is AngleUnit.MM_OVER_METER:
        return np.arctan(normalize_numeric_type(value) / 1000.0)

    # ----------------------------------------------------------------------------------
    if unit is AngleUnit.DMS:
        array = as_array(value)

        if array.ndim not in (1, 2) or array.shape[-1] != 3:
            raise ValueError(
                "DMS input must have shape (3,) or (N, 3): "
                "[degrees, arcminutes, arcseconds]."
            )

        if np.any(array[..., 1:] < 0.0):
            raise ValueError("DMS minutes and seconds must be non-negative.")

        sign = np.where(array[..., 0] < 0.0, -1.0, 1.0)
        degrees = np.abs(array[..., 0]) + array[..., 1] / 60.0 + array[..., 2] / 3600.0

        return sign * degrees * ANGLE_FACTORS[AngleUnit.DEGREE]

    # ----------------------------------------------------------------------------------
    return normalize_numeric_type(value) * ANGLE_FACTORS[unit]


# ======================================================================================
def _radians_to_angle(radians: object, unit: AngleUnit) -> float | FloatArray:
    """Convert an angle in radians to the specified target unit."""

    if unit is AngleUnit.MM_OVER_METER:
        return 1000.0 * np.tan(normalize_numeric_type(radians))

    # ----------------------------------------------------------------------------------
    if unit is AngleUnit.DMS:
        degrees = as_array(radians) / ANGLE_FACTORS[AngleUnit.DEGREE]

        sign = np.where(degrees < 0.0, -1.0, 1.0)
        absolute_degrees = np.abs(degrees)

        whole_degrees = np.floor(absolute_degrees)
        minutes_total = (absolute_degrees - whole_degrees) * 60.0
        whole_minutes = np.floor(minutes_total)
        seconds = (minutes_total - whole_minutes) * 60.0

        whole_degrees = sign * whole_degrees

        return np.stack([whole_degrees, whole_minutes, seconds], axis=-1)

    # ----------------------------------------------------------------------------------
    return normalize_numeric_type(radians) / ANGLE_FACTORS[unit]


# ======================================================================================
def _convert_angle(
    value: object, source: AngleUnit, target: AngleUnit
) -> float | FloatArray:
    """Convert an angle between supported units."""
    radians = _angle_to_radians(value, source)
    return _radians_to_angle(radians, target)


# ======================================================================================
def _convert_distance(
    value: object, source: DistanceUnit, target: DistanceUnit
) -> float | FloatArray:
    """Convert a distance between supported units."""
    factor = DISTANCE_FACTORS[source] / DISTANCE_FACTORS[target]

    return normalize_numeric_type(value) * factor


# ======================================================================================
# Internal conversion interface
# ======================================================================================
class _Conversion:
    """Callable conversion between a source unit and a target unit."""

    def __init__(
        self, source: AngleUnit | DistanceUnit, target: AngleUnit | DistanceUnit
    ) -> None:

        if type(source) is not type(target):
            raise ValueError(
                f"Cannot convert {source!r} to {target!r}: "
                "units belong to different dimensions."
            )

        # ------------------------------------------------------------------------------
        self._source = source
        self._target = target

    # ==================================================================================
    @overload
    def __call__(self, value: Number | np.number) -> float: ...

    @overload
    def __call__(self, value: ArrayLike) -> FloatArray: ...

    def __call__(self, value: object) -> float | FloatArray:
        """Convert value from the source unit to the target unit."""
        if isinstance(self._source, AngleUnit) and isinstance(self._target, AngleUnit):
            return _convert_angle(value, self._source, self._target)
        # ------------------------------------------------------------------------------
        if isinstance(self._source, DistanceUnit) and isinstance(
            self._target, DistanceUnit
        ):
            return _convert_distance(value, self._source, self._target)
        # ------------------------------------------------------------------------------
        raise ValueError(
            f"Cannot convert {self._source!r} to {self._target!r}: "
            "units belong to different dimensions."
        )


# ======================================================================================
# Angle target units
# ======================================================================================
class _AngleTo:
    """Expose supported target units for angle conversions."""

    def __init__(self, source: AngleUnit) -> None:
        self._source = source

    def __getattr__(self, name: str) -> _Conversion:
        unit = _resolve_unit(name)

        if isinstance(unit, DistanceUnit):
            raise AttributeError(f"Cannot convert angle to length unit: {name!r}.")

        raise AttributeError(f"Unknown angular unit: {name!r}.")

    # ==================================================================================
    # Properties for IDE autocomplete
    # ==================================================================================
    @property
    def radian(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.RADIAN)

    @property
    def milliradian(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.MILLIRADIAN)

    @property
    def mm_over_meter(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.MM_OVER_METER)

    # ==================================================================================
    @property
    def gon(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.GON)

    @property
    def grad(self) -> _Conversion:
        """Alias property for gon."""
        return _Conversion(self._source, AngleUnit.GON)

    @property
    def c(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.C)

    @property
    def cc(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.CC)

    # ==================================================================================
    @property
    def degree(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.DEGREE)

    @property
    def arcminute(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.ARCMINUTE)

    @property
    def arcsecond(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.ARCSECOND)

    @property
    def dms(self) -> _Conversion:
        return _Conversion(self._source, AngleUnit.DMS)


# ======================================================================================
# Length target units
# ======================================================================================
class _LengthTo:
    """Expose supported target units for length conversions."""

    def __init__(self, source: DistanceUnit) -> None:
        self._source = source

    def __getattr__(self, name: str) -> _Conversion:
        unit = _resolve_unit(name)

        if isinstance(unit, AngleUnit):
            raise AttributeError(f"Cannot convert length to angular unit: {name!r}.")

        raise AttributeError(f"Unknown length unit: {name!r}.")

    # ==================================================================================
    # Properties for IDE autocomplete
    # ==================================================================================
    @property
    def kilometer(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.KILOMETER)

    @property
    def hectometer(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.HECTOMETER)

    @property
    def meter(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.METER)

    @property
    def decimeter(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.DECIMETER)

    @property
    def centimeter(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.CENTIMETER)

    @property
    def millimeter(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.MILLIMETER)

    @property
    def micrometer(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.MICROMETER)

    # ==================================================================================
    @property
    def mile(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.MILE)

    @property
    def yard(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.YARD)

    @property
    def feet(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.FEET)

    @property
    def inch(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.INCH)

    @property
    def thou(self) -> _Conversion:
        return _Conversion(self._source, DistanceUnit.THOU)


# ======================================================================================
# Angle source units
# ======================================================================================
class _AngleFrom:
    """Expose target-unit conversions for a source angle unit."""

    def __init__(self, source: AngleUnit) -> None:
        self._source = source

    # ==================================================================================
    @property
    def to(self) -> _AngleTo:
        return _AngleTo(self._source)


# ======================================================================================
# Length source units
# ======================================================================================
class _LengthFrom:
    """Expose target-unit conversions for a source length unit."""

    def __init__(self, source: DistanceUnit) -> None:
        self._source = source

    # ==================================================================================
    @property
    def to(self) -> _LengthTo:
        return _LengthTo(self._source)


# ======================================================================================
# Public conversion entry point
# ======================================================================================
class _From:
    """Expose all supported source units for unit conversion."""

    def __getattr__(self, name: str) -> _AngleFrom | _LengthFrom:
        unit = _resolve_unit(name)

        if isinstance(unit, AngleUnit):
            return _AngleFrom(unit)

        if isinstance(unit, DistanceUnit):
            return _LengthFrom(unit)

        raise AttributeError(f"Unknown unit: {name!r}.")

    # ==================================================================================
    # Angle sources
    # ==================================================================================
    @property
    def radian(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.RADIAN)

    @property
    def milliradian(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.MILLIRADIAN)

    @property
    def mm_over_meter(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.MM_OVER_METER)

    # ==================================================================================
    @property
    def gon(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.GON)

    @property
    def grad(self) -> _AngleFrom:
        """Alias property for gon."""
        return _AngleFrom(AngleUnit.GON)

    @property
    def c(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.C)

    @property
    def cc(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.CC)

    # ==================================================================================
    @property
    def degree(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.DEGREE)

    @property
    def arcminute(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.ARCMINUTE)

    @property
    def arcsecond(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.ARCSECOND)

    @property
    def dms(self) -> _AngleFrom:
        return _AngleFrom(AngleUnit.DMS)

    # ==================================================================================
    # Length sources
    # ==================================================================================
    @property
    def kilometer(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.KILOMETER)

    @property
    def hectometer(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.HECTOMETER)

    @property
    def meter(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.METER)

    @property
    def decimeter(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.DECIMETER)

    @property
    def centimeter(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.CENTIMETER)

    @property
    def millimeter(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.MILLIMETER)

    @property
    def micrometer(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.MICROMETER)

    # ==================================================================================
    @property
    def mile(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.MILE)

    @property
    def yard(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.YARD)

    @property
    def feet(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.FEET)

    @property
    def inch(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.INCH)

    @property
    def thou(self) -> _LengthFrom:
        return _LengthFrom(DistanceUnit.THOU)


# ======================================================================================
# Public API
# ======================================================================================
convert = _From()
