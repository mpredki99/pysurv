# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.utils.aliased_str_enum import AliasedStrEnum


# ======================================================================================
# Distance units enumeration
# ======================================================================================
class DistanceUnit(AliasedStrEnum):
    """ "PySurveyor distance units enumeration class."""

    __to_lower__ = True

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
class AngleUnit(AliasedStrEnum):
    """ "PySurveyor angle units enumeration class."""

    __to_lower__ = True

    # ==================================================================================
    # Angle units
    # ==================================================================================
    RADIAN = "radian", "radians", "rad", "rads"
    MILIRADIAN = "milliradian", "milliradians", "mrad", "mrads"
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
