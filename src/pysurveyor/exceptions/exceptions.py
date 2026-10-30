# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.


# ======================================================================================
class PySurveyorError(Exception):
    """Base class for all PySurveyor exceptions."""

    pass


# ======================================================================================
class ValidationError(PySurveyorError):
    """Raised by a PySurveyor validator when validation fails."""

    pass


# ======================================================================================
class ParseError(PySurveyorError):
    """Raised when a parsing error occurs."""

    pass
