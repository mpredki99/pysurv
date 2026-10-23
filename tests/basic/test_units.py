"""Tests for the pysurveyo units."""

import pytest

from pysurveyor.basic import AngleUnit, DistanceUnit


# ======================================================================================
# Tests for AngleUnit and DistanceUnit
# ======================================================================================
@pytest.mark.parametrize(
    ("enum_class", "canonical", "aliases"),
    [
        # Metric distance units --------------------------------------------------------
        (DistanceUnit, "kilometer", ("kilometers", "km")),
        (DistanceUnit, "hectometer", ("hectometers", "hm")),
        (DistanceUnit, "meter", ("meters", "m")),
        (DistanceUnit, "decimeter", ("decimeters", "dm")),
        (DistanceUnit, "centimeter", ("centimeters", "cm")),
        (DistanceUnit, "millimeter", ("millimeters", "mm")),
        (DistanceUnit, "micrometer", ("micrometers", "um")),
        # Imperial distance units ------------------------------------------------------
        (DistanceUnit, "mile", ("miles", "mi")),
        (DistanceUnit, "yard", ("yards", "yd")),
        (DistanceUnit, "feet", ("foot", "ft")),
        (DistanceUnit, "inch", ("inches", "in")),
        (DistanceUnit, "thou", ("thous", "th")),
        # Angle units ------------------------------------------------------------------
        (AngleUnit, "radian", ("radians", "rad", "rads")),
        (AngleUnit, "milliradian", ("milliradians", "mrad", "mrads")),
        (AngleUnit, "gon", ("gons", "grad", "grads", "g")),
        (AngleUnit, "c", ()),
        (AngleUnit, "cc", ()),
        (AngleUnit, "degree", ("degrees", "deg", "degs")),
        (AngleUnit, "arcminute", ("arcminutes", "arcmin", "arcmins", "min", "mins")),
        (AngleUnit, "arcsecond", ("arcseconds", "arcsec", "arcsecs", "sec", "secs")),
        (
            AngleUnit,
            "mm_over_meter",
            (
                "mm_over_meter",
                "mm_over_m",
                "mm/m",
                "millimeter per meter",
                "millimeters per meter",
            ),
        ),
        (
            AngleUnit,
            "dms",
            (
                "degree minute second",
                "degrees minutes seconds",
                "deg arcmin arcsec",
                "deg min sec",
            ),
        ),
    ],
)
def test_canonical_values_and_aliases(enum_class, canonical, aliases):
    member = enum_class(canonical)

    assert member.value == canonical
    assert str(member) == canonical

    for alias in aliases:
        assert enum_class(alias) is member


# ======================================================================================
@pytest.mark.parametrize(
    ("enum_class", "canonical", "aliases"),
    [
        # Metric distance units --------------------------------------------------------
        (DistanceUnit, "kilometer", ("kilometers", "km")),
        (DistanceUnit, "hectometer", ("hectometers", "hm")),
        (DistanceUnit, "meter", ("meters", "m")),
        (DistanceUnit, "decimeter", ("decimeters", "dm")),
        (DistanceUnit, "centimeter", ("centimeters", "cm")),
        (DistanceUnit, "millimeter", ("millimeters", "mm")),
        (DistanceUnit, "micrometer", ("micrometers", "um")),
        # Imperial distance units ------------------------------------------------------
        (DistanceUnit, "mile", ("miles", "mi")),
        (DistanceUnit, "yard", ("yards", "yd")),
        (DistanceUnit, "feet", ("foot", "ft")),
        (DistanceUnit, "inch", ("inches", "in")),
        (DistanceUnit, "thou", ("thous", "th")),
        # Angle units ------------------------------------------------------------------
        (AngleUnit, "radian", ("radians", "rad", "rads")),
        (AngleUnit, "milliradian", ("milliradians", "mrad", "mrads")),
        (AngleUnit, "gon", ("gons", "grad", "grads", "g")),
        (AngleUnit, "c", ()),
        (AngleUnit, "cc", ()),
        (AngleUnit, "degree", ("degrees", "deg", "degs")),
        (AngleUnit, "arcminute", ("arcminutes", "arcmin", "arcmins", "min", "mins")),
        (AngleUnit, "arcsecond", ("arcseconds", "arcsec", "arcsecs", "sec", "secs")),
        (
            AngleUnit,
            "mm_over_meter",
            (
                "mm_over_meter",
                "mm_over_m",
                "mm/m",
                "millimeter per meter",
                "millimeters per meter",
            ),
        ),
        (
            AngleUnit,
            "dms",
            (
                "degree minute second",
                "degrees minutes seconds",
                "deg arcmin arcsec",
                "deg min sec",
            ),
        ),
    ],
)
def test_case_insensitive_lookup(enum_class, canonical, aliases):
    member = enum_class(canonical)

    for value in (canonical, *aliases):
        assert enum_class(value.upper()) is member
        assert enum_class(value.capitalize()) is member


# ======================================================================================
@pytest.mark.parametrize(
    ("enum_class", "invalid_value"),
    [
        (DistanceUnit, "mile-per-hour"),
        (DistanceUnit, "unknown"),
        (AngleUnit, "unknown"),
        (AngleUnit, "degreee"),
        (DistanceUnit, ""),
        (AngleUnit, ""),
        (DistanceUnit, None),
        (AngleUnit, None),
        (DistanceUnit, 123),
        (AngleUnit, 123),
    ],
)
def test_invalid_values(enum_class, invalid_value):
    with pytest.raises(ValueError):
        enum_class(invalid_value)
