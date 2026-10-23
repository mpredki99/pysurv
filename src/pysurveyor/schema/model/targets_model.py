# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.basic.units import AngleUnit, DistanceUnit
from pysurveyor.validators.comparison_validators import GreaterOrEqual
from pysurveyor.validators.conversion_validators import IsNumeric, IsText
from pysurveyor.validators.utils.constants import COMMENT


def get_targets_model_data() -> dict:
    """Returns a dictionary containing the schema for the targets model."""

    return {
        "trg_id": [True, None, IsText(ignore=[]), "Target identifier"],
        # ------------------------------------------------------------------------------
        "trg_cst": [False, DistanceUnit.METER, IsNumeric(), "Target constant"],
        "u_trg_cst": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in target constant",
        ],
        # ------------------------------------------------------------------------------
        "u_trg_ctr": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in target centering",
        ],
        # ------------------------------------------------------------------------------
        "u_trg_tilt": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in target tilt",
        ],
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
    }
