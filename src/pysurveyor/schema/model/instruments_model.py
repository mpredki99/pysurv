# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.basic.units import AngleUnit, DistanceUnit
from pysurveyor.validators.comparison_validators import GreaterOrEqual
from pysurveyor.validators.conversion_validators import IsNumeric, IsText
from pysurveyor.validators.utils.constants import COMMENT


def get_instruments_model_data() -> dict:
    """Returns a dictionary containing the schema for the instruments model."""

    return {
        "instr_id": [True, None, IsText(ignore=[]), "Instrument identifier"],
        # ------------------------------------------------------------------------------
        "instr_b": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Instrument trunnion error angle",
        ],
        "u_instr_b": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in instrument trunnion error angle",
        ],
        # ------------------------------------------------------------------------------
        "instr_c": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Instrument collimation error angle",
        ],
        "u_instr_c": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in instrument collimation error angle",
        ],
        # ------------------------------------------------------------------------------
        "instr_i": [
            False,
            AngleUnit.RADIAN,
            IsNumeric(),
            "Instrument index error angle",
        ],
        "u_instr_i": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in instrument index error angle",
        ],
        # ------------------------------------------------------------------------------
        "u_instr_hz": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in horizontal direction reading",
        ],
        # ------------------------------------------------------------------------------
        "u_instr_v": [
            False,
            AngleUnit.RADIAN,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in vertical angle reading",
        ],
        # ------------------------------------------------------------------------------
        "u_instr_sd": [
            False,
            DistanceUnit.METER,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in spatial distance measurement",
        ],
        "u_instr_sd_ppm": [
            False,
            None,
            IsNumeric() & GreaterOrEqual(0, ignore=[None, COMMENT, -1]),
            "Uncertainty in spatial distance measurement unit per million",
        ],
    }
