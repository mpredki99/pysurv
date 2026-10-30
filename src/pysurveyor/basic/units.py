# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.utils.aliased_str_enum import AliasedStrEnum


# ======================================================================================
# Base class for ``PySurveyor`` units enumeration
# ======================================================================================
class Unit(AliasedStrEnum):
    """Base class for ``PySurveyor`` units enumeration classes."""

    __to_lower__ = True


# ======================================================================================
# Distance units enumeration
# ======================================================================================
class DistanceUnit(Unit):
    """``PySurveyor`` distance units enumeration class."""

    # ==================================================================================
    # Metric distance units
    # ==================================================================================
    KILOMETER = "kilometer", "kilometers", "km"
    HECTOMETER = "hectometer", "hectometers", "hm"
    METER = "meter", "meters", "m"
    DECIMETER = "decimeter", "decimeters", "dm"
    CENTIMETER = "centimeter", "centimeters", "cm"
    MILLIMETER = "millimeter", "millimeters", "mm"
    MICROMETER = "micrometer", "micrometers", "um"
    # ==================================================================================
    # Imperial distance units
    # ==================================================================================
    MILE = "mile", "miles", "mi"
    YARD = "yard", "yards", "yd"
    FEET = "feet", "foot", "ft"
    INCH = "inch", "inches", "in"
    THOU = "thou", "thous", "th"


# ======================================================================================
# Angle units enumeration
# ======================================================================================
class AngleUnit(Unit):
    """``PySurveyor`` angle units enumeration class."""

    # ==================================================================================
    # Angle units
    # ==================================================================================
    RADIAN = "radian", "radians", "rad", "rads"
    MILLIRADIAN = "milliradian", "milliradians", "mrad", "mrads"
    GON = "gon", "gons", "grad", "grads", "g"
    C = "c"
    CC = "cc"
    DEGREE = "degree", "degrees", "deg", "degs"
    ARCMINUTE = "arcminute", "arcminutes", "arcmin", "arcmins", "min", "mins"
    ARCSECOND = "arcsecond", "arcseconds", "arcsec", "arcsecs", "sec", "secs"
    MM_OVER_METER = (
        "mm_over_meter",
        "mm_over_m",
        "mm/m",
        "millimeter per meter",
        "millimeters per meter",
    )
    DMS = (
        "dms",
        "degree minute second",
        "degrees minutes seconds",
        "deg arcmin arcsec",
        "deg min sec",
    )


# ======================================================================================
# Parse unit
# ======================================================================================
def parse_unit(unit: str) -> Unit:
    """Parse a unit string and return the corresponding ``Unit`` enumeration member.

    Parameters
    ----------
    unit : str
        The unit string to parse.

    Returns
    -------
    Unit
        The corresponding ``Unit`` enumeration member.
    """
    try:
        return DistanceUnit(unit)
    except ValueError:
        pass

    try:
        return AngleUnit(unit)
    except ValueError:
        pass

    raise ValueError(f"Unknown unit: {unit!r}")
