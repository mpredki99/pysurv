# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.basic.units import AngleUnit, DistanceUnit
from pysurveyor.validators.comparison_validators import GreaterOrEqual
from pysurveyor.validators.conversion_validators import IsNumeric, IsText
from pysurveyor.validators.utils.constants import COMMENT


def get_stations_model_data() -> dict:
    """Returns a dictionary containing the schema for the stations model."""

    return {
        "stn_id": [True, None, IsText(ignore=[]), "Station identifier"],
        # ------------------------------------------------------------------------------
        "stn_pt": [True, None, IsText(ignore=[]), "Station base point identifier"],
        # ------------------------------------------------------------------------------
        "instr_id": [True, None, IsText(ignore=[]), "Instrument identifier"],
        # ------------------------------------------------------------------------------
        "stn_h": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "Station height above the base point",
        ],
        "u_stn_h": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in station height",
        ],
        # ------------------------------------------------------------------------------
        "stn_rx": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Instrument rotation about the local X-axis",
        ],
        "u_stn_rx": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in instrument rotation about the local X-axis",
        ],
        # ------------------------------------------------------------------------------
        "stn_ry": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Instrument rotation about the local Y-axis",
        ],
        "u_stn_ry": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in instrument rotation about the local Y-axis",
        ],
        # ------------------------------------------------------------------------------
        "stn_rz": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Instrument orientation with respect to GLOBAL y-axis",
        ],
        "u_stn_rz": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in instrument orientation with respect to GLOBAL y-axis",
        ],
    }
