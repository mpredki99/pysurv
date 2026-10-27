# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.basic.units import AngleUnit, DistanceUnit
from pysurveyor.validators.comparison_validators import GreaterOrEqual
from pysurveyor.validators.conversion_validators import IsNumeric, IsText
from pysurveyor.validators.utils.constants import COMMENT


def get_control_points_model_data() -> dict:
    """Returns a dictionary containing the schema for the controls model."""

    return {
        "id": [True, None, IsText(ignore=[]), "Control point identifier"],
        # ------------------------------------------------------------------------------
        "x": [False, DistanceUnit.METER, IsNumeric(), "X coordinate of the point"],
        "u_x": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0),
            "Uncertainty in X coordinate",
        ],
        # ------------------------------------------------------------------------------
        "y": [False, DistanceUnit.METER, IsNumeric(), "Y coordinate of the point"],
        "u_y": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in Y coordinate",
        ],
        # ------------------------------------------------------------------------------
        "z": [False, DistanceUnit.METER, IsNumeric(), "Z coordinate of the point"],
        "u_z": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in Z coordinate",
        ],
        # ------------------------------------------------------------------------------
        "xi": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "South-North component of the deflection of the vertical",
        ],
        "u_xi": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in xi",
        ],
        # ------------------------------------------------------------------------------
        "eta": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "West-East component of the deflection of the vertical",
        ],
        "u_eta": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in eta",
        ],
        # ------------------------------------------------------------------------------
        "n": [
            False,
            DistanceUnit.METER,
            IsNumeric(),
            "Undulation or height anomaly at the point",
        ],
        "u_n": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in n",
        ],
    }
