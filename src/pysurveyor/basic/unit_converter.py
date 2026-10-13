"""Unit conversion utilities for angles and lengths.

This module provides an API for converting between angular and linear 
units using NumPy-compatible numerical operations.

Supported angle units
---------------------
- `rad`: radians
- `mrad`: milliradians
- `gon` and `grad`: gradians (400 gon per full revolution)
- `c`: centesimal angular minutes (1/100 gon)
- `cc`: centesimal angular seconds (1/10000 gon)
- `degrees`: decimal degrees
- `min`: angular minutes (1/60 degree)
- `sec`: angular seconds (1/3600 degree)
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
ValueError
    If a source or target unit is unsupported, or if DMS input has an
    invalid shape or negative minutes or seconds.
"""

from __future__ import annotations

from typing import overload

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


# ======================================================================================
# Unit conversion factors
# ======================================================================================
ANGLE_FACTORS: dict[str, float] = {
    "rad":     1.0,
    "mrad":    0.001,
    "gon":     np.pi / 200.0,
    "c":       np.pi / 20000.0,
    "cc":      np.pi / 2000000.0,
    "degrees": np.pi / 180.0,
    "min":     np.pi / 10800.0,
    "sec":     np.pi / 648000.0,
}


# ======================================================================================
LENGTH_FACTORS: dict[str, float] = {
    "kilometer":   1000.0,
    "hectometer":  100.0,
    "meter":       1.0,
    "decimeter":   0.1,
    "centimeter":  0.01,
    "millimeters": 0.001,
    "micrometer":  0.000001,
    # ------------------------------
    "mile": 1609.344,
    "yard": 0.9144,
    "feet": 0.3048,
    "inch": 0.0254,
    "thou": 0.0000254,
}


# ======================================================================================
# Conversion helper functions
# ======================================================================================

def _is_scalar(value: object) -> bool:
    """Return whether value is a supported scalar numeric type."""
    return isinstance(value, (int, float, np.number)) or (
        isinstance(value, np.ndarray) and value.ndim == 0
    )


# ======================================================================================
def _as_array(value: object) -> FloatArray:
    """Convert value to a NumPy array with float64 precision."""
    return np.asarray(value, dtype=np.float64)


# ======================================================================================
def _angle_to_radians(value: object, unit: str) -> float | FloatArray:
    """Convert an angle from the specified unit to radians."""

    if unit == "mm_over_meter":
        if _is_scalar(value):
            return float(np.arctan(float(value) / 1000.0))
        return np.arctan(_as_array(value) / 1000.0)

    # ----------------------------------------------------------------------------------
    if unit == "dms":
        array = _as_array(value)

        if array.ndim not in (1, 2) or array.shape[-1] != 3:
            raise ValueError(
                "DMS input must have shape (3,) or (N, 3): "
                "[degrees, minutes, seconds]."
            )

        if np.any(array[..., 1:] < 0.0):
            raise ValueError("DMS minutes and seconds must be non-negative.")

        sign = np.where(array[..., 0] < 0.0, -1.0, 1.0)
        degrees = (
            np.abs(array[..., 0])
            + array[..., 1] / 60.0
            + array[..., 2] / 3600.0
        )

        return sign * degrees * ANGLE_FACTORS["degrees"]

    # ----------------------------------------------------------------------------------
    factor = ANGLE_FACTORS[unit]

    if _is_scalar(value):
        return float(value) * factor

    return _as_array(value) * factor


# ======================================================================================
def _radians_to_angle(radians: float | FloatArray, unit: str) -> float | FloatArray:
    """Convert an angle in radians to the specified target unit."""

    if unit == "mm_over_meter":
        if isinstance(radians, float):
            return float(1000.0 * np.tan(radians))

        return 1000.0 * np.tan(np.asarray(radians, dtype=np.float64))

    # ----------------------------------------------------------------------------------
    if unit == "dms":
        degrees = np.asarray(radians, dtype=np.float64) / ANGLE_FACTORS["degrees"]

        sign = np.where(degrees < 0.0, -1.0, 1.0)
        absolute_degrees = np.abs(degrees)

        whole_degrees = np.floor(absolute_degrees)
        minutes_total = (absolute_degrees - whole_degrees) * 60.0
        whole_minutes = np.floor(minutes_total)
        seconds = (minutes_total - whole_minutes) * 60.0

        whole_degrees = sign * whole_degrees

        return np.stack([whole_degrees, whole_minutes, seconds], axis=-1)

    # ----------------------------------------------------------------------------------
    factor = ANGLE_FACTORS[unit]

    if _is_scalar(radians):
        return float(radians) / factor

    return np.asarray(radians, dtype=np.float64) / factor


# ======================================================================================
def _convert_angle(value: object, source: str, target: str) -> float | FloatArray:
    """Convert an angle between supported units."""
    radians = _angle_to_radians(value, source)
    return _radians_to_angle(radians, target)


# ======================================================================================
def _convert_length(value: object, source: str, target: str) -> float | FloatArray:
    """Convert a length between supported units."""
    factor = LENGTH_FACTORS[source] / LENGTH_FACTORS[target]

    if _is_scalar(value):
        return float(value) * factor

    return _as_array(value) * factor


# ======================================================================================
# Internal conversion interface
# ======================================================================================
class _Conversion:
    """Callable conversion between a source unit and a target unit."""

    def __init__(self, source: str, target: str, dimension: str) -> None:
        angle_units = set(ANGLE_FACTORS) | {"dms", "mm_over_meter"}
        length_units = set(LENGTH_FACTORS)

        if dimension == "angle":
            valid_units = angle_units
        elif dimension == "length":
            valid_units = length_units
        else:
            raise ValueError(f"Unsupported dimension: {dimension!r}")

        if source not in valid_units:
            raise ValueError(f"Source unit {source!r} is not a valid {dimension} unit.")

        if target not in valid_units:
            raise ValueError(f"Target unit {target!r} is not a valid {dimension} unit.")

        # ------------------------------------------------------------------------------
        self._source = source
        self._target = target
        self._dimension = dimension

    # ==================================================================================
    @overload
    def __call__(self, value: int | float | np.number) -> float: ...

    @overload
    def __call__(self, value: ArrayLike) -> FloatArray: ...

    def __call__(self, value: object) -> float | FloatArray:
        """Convert value from the source unit to the target unit."""
        if self._dimension == "angle":
            return _convert_angle(value, self._source, self._target)

        return _convert_length(value, self._source, self._target)


# ======================================================================================
# Angle target units
# ======================================================================================
class _AngleTo:
    """Expose supported target units for angle conversions."""

    def __init__(self, source: str) -> None:
        self._source = source

    def __getattr__(self, name: str) -> _Conversion:
        if name in LENGTH_FACTORS:
            raise AttributeError(f"Cannot convert angle to length unit: {name!r}.")

        raise AttributeError(f"Invalid angular unit: {name!r}.")

    # ==================================================================================
    # Properties for IDE autocomplete
    # ==================================================================================
    @property
    def rad(self) -> _Conversion:
        return _Conversion(self._source, "rad", "angle")

    @property
    def mrad(self) -> _Conversion:
        return _Conversion(self._source, "mrad", "angle")

    @property
    def mm_over_meter(self) -> _Conversion:
        return _Conversion(self._source, "mm_over_meter", "angle")

    # ==================================================================================
    @property
    def gon(self) -> _Conversion:
        return _Conversion(self._source, "gon", "angle")

    @property
    def grad(self) -> _Conversion:
        """Alias property for gon."""
        return _Conversion(self._source, "gon", "angle")

    @property
    def c(self) -> _Conversion:
        return _Conversion(self._source, "c", "angle")

    @property
    def cc(self) -> _Conversion:
        return _Conversion(self._source, "cc", "angle")

    # ==================================================================================
    @property
    def degrees(self) -> _Conversion:
        return _Conversion(self._source, "degrees", "angle")

    @property
    def min(self) -> _Conversion:
        return _Conversion(self._source, "min", "angle")

    @property
    def sec(self) -> _Conversion:
        return _Conversion(self._source, "sec", "angle")

    @property
    def dms(self) -> _Conversion:
        return _Conversion(self._source, "dms", "angle")


# ======================================================================================
# Length target units
# ======================================================================================
class _LengthTo:
    """Expose supported target units for length conversions."""

    def __init__(self, source: str) -> None:
        self._source = source

    def __getattr__(self, name: str) -> _Conversion:
        if name in set(ANGLE_FACTORS) | {"dms", "mm_over_meter", "grad"}:
            raise AttributeError(f"Cannot convert length to angular unit: {name!r}.")

        raise AttributeError(f"Invalid length unit: {name!r}.")

    # ==================================================================================
    # Properties for IDE autocomplete
    # ==================================================================================
    @property
    def kilometer(self) -> _Conversion:
        return _Conversion(self._source, "kilometer", "length")

    @property
    def hectometer(self) -> _Conversion:
        return _Conversion(self._source, "hectometer", "length")
    
    @property
    def meter(self) -> _Conversion:
        return _Conversion(self._source, "meter", "length")

    @property
    def decimeter(self) -> _Conversion:
        return _Conversion(self._source, "decimeter", "length")

    @property
    def centimeter(self) -> _Conversion:
        return _Conversion(self._source, "centimeter", "length")

    @property
    def millimeters(self) -> _Conversion:
        return _Conversion(self._source, "millimeters", "length")

    @property
    def micrometer(self) -> _Conversion:
        return _Conversion(self._source, "micrometer", "length")

    # ==================================================================================
    @property
    def mile(self) -> _Conversion:
        return _Conversion(self._source, "mile", "length")

    @property
    def yard(self) -> _Conversion:
        return _Conversion(self._source, "yard", "length")

    @property
    def feet(self) -> _Conversion:
        return _Conversion(self._source, "feet", "length")
    
    @property
    def inch(self) -> _Conversion:
        return _Conversion(self._source, "inch", "length")

    @property
    def thou(self) -> _Conversion:
        return _Conversion(self._source, "thou", "length")


# ======================================================================================
# Angle source units
# ======================================================================================
class _AngleFrom:
    """Expose target-unit conversions for a source angle unit."""

    def __init__(self, source: str) -> None:
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

    def __init__(self, source: str) -> None:
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
        if name in set(ANGLE_FACTORS) | {"dms", "mm_over_meter", "grad"}:
            return _AngleFrom("gon" if name == "grad" else name)

        if name in LENGTH_FACTORS:
            return _LengthFrom(name)

        raise AttributeError(f"Invalid unit: {name!r}.")

    # ==================================================================================
    # Angle sources
    # ==================================================================================
    @property
    def rad(self) -> _AngleFrom:
        return _AngleFrom("rad")

    @property
    def mrad(self) -> _AngleFrom:
        return _AngleFrom("mrad")

    @property
    def mm_over_meter(self) -> _AngleFrom:
        return _AngleFrom("mm_over_meter")

    # ==================================================================================
    @property
    def gon(self) -> _AngleFrom:
        return _AngleFrom("gon")

    @property
    def grad(self) -> _AngleFrom:
        """Alias property for gon."""
        return _AngleFrom("gon")

    @property
    def c(self) -> _AngleFrom:
        return _AngleFrom("c")

    @property
    def cc(self) -> _AngleFrom:
        return _AngleFrom("cc")

    # ==================================================================================
    @property
    def degrees(self) -> _AngleFrom:
        return _AngleFrom("degrees")

    @property
    def min(self) -> _AngleFrom:
        return _AngleFrom("min")

    @property
    def sec(self) -> _AngleFrom:
        return _AngleFrom("sec")

    @property
    def dms(self) -> _AngleFrom:
        return _AngleFrom("dms")

    # ==================================================================================
    # Length sources
    # ==================================================================================
    @property
    def kilometer(self) -> _LengthFrom:
        return _LengthFrom("kilometer")

    @property
    def hectometer(self) -> _LengthFrom:
        return _LengthFrom("hectometer")

    @property
    def meter(self) -> _LengthFrom:
        return _LengthFrom("meter")

    @property
    def decimeter(self) -> _LengthFrom:
        return _LengthFrom("decimeter")

    @property
    def centimeter(self) -> _LengthFrom:
        return _LengthFrom("centimeter")

    @property
    def millimeters(self) -> _LengthFrom:
        return _LengthFrom("millimeters")

    @property
    def micrometer(self) -> _LengthFrom:
        return _LengthFrom("micrometer")

    # ==================================================================================
    @property
    def mile(self) -> _LengthFrom:
        return _LengthFrom("mile")

    @property
    def yard(self) -> _LengthFrom:
        return _LengthFrom("yard")

    @property
    def feet(self) -> _LengthFrom:
        return _LengthFrom("feet")
    
    @property
    def inch(self) -> _LengthFrom:
        return _LengthFrom("inch")

    @property
    def thou(self) -> _LengthFrom:
        return _LengthFrom("thou")


# ======================================================================================
# Public API
# ======================================================================================
convert = _From()
