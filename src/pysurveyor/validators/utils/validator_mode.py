# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from enum import StrEnum


# ======================================================================================
# Validator modes supported by PySurveyor validators
# ======================================================================================
class ValidatorMode(StrEnum):
    """Define how a validator handles validation failures.

    Members
    -------
    RETURN
        Return a boolean mask indicating valid and invalid values.
    RAISE
        Raise a validation error when invalid values are found.
    """

    RETURN = "return"
    RAISE = "raise"
