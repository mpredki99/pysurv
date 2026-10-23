# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.basic.units import AngleUnit, DistanceUnit
from pysurveyor.validators.comparison_validators import GreaterOrEqual
from pysurveyor.validators.conversion_validators import IsNumeric, IsText
from pysurveyor.validators.utils.constants import COMMENT


def get_measurements_model_data() -> dict:
    """Returns a dictionary containing the schema for the measurements model."""

    return {
        "stn_id": [True, None, IsText(ignore=[]), "Station identifier"],
        # ------------------------------------------------------------------------------
        "trg_pt": [True, None, IsText(ignore=[]), "Target point identifier"],
        # ------------------------------------------------------------------------------
        "trg_id": [True, None, IsText(ignore=[]), "Target identifier"],
        # ------------------------------------------------------------------------------
        "trg_h": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "Target height above the base point",
        ],
        "u_trg_h": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in target height",
        ],
        "u_trg_ctr": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in target centering",
        ],
        "u_trg_tilt": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in target tilt",
        ],
        # ------------------------------------------------------------------------------
        "sd": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "Slope distance from the station to the target",
        ],
        "u_sd": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in slope distance",
        ],
        # ------------------------------------------------------------------------------
        "hd": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "Horizontal distance from the station to the target",
        ],
        "u_hd": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in horizontal distance",
        ],
        # ------------------------------------------------------------------------------
        "vd": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "Vertical distance from the station to the target",
        ],
        "u_vd": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in vertical distance",
        ],
        # ------------------------------------------------------------------------------
        "dx": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "X-axis component of the GNSS vector from the station to the target",
        ],
        "u_dx": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in X-axis component of the GNSS vector",
        ],
        # ------------------------------------------------------------------------------
        "dy": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "Y-axis component of the GNSS vector from the station to the target",
        ],
        "u_dy": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in Y-axis component of the GNSS vector",
        ],
        # ------------------------------------------------------------------------------
        "dz": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "Z-axis component of the GNSS vector from the station to the target",
        ],
        "u_dz": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in Z-axis component of the GNSS vector",
        ],
        # ------------------------------------------------------------------------------
        "az": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Azimuth from the station to the target",
        ],
        "u_az": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in azimuth",
        ],
        # ------------------------------------------------------------------------------
        "hz": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Horizontal direction reading from the station to the target",
        ],
        "u_hz": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in horizontal direction reading",
        ],
        # ------------------------------------------------------------------------------
        "vz": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Vertical zenithal direction reading from the station to the target",
        ],
        "u_vz": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in vertical zenithal direction reading",
        ],
        # ------------------------------------------------------------------------------
        "vh": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Vertical direction reading from horizontal plane to the target",
        ],
        "u_vh": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in vertical direction reading from horizontal plane",
        ],
    }
