# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from enum import StrEnum


# ======================================================================================
# Validation modes supported by PySurveyor schemas
# ======================================================================================
class ValidationMode(StrEnum):
    """Define how a schema validates its associated data table.

    The enum inherits from :class:`StrEnum`, so each member can be used
    wherever the corresponding lowercase string value is expected.

    Members
    -------
    DISABLED
        Skip validation.
    EAGER
        Raise a validation error as soon as the first failure is found.
    LAZY
        Collect validation failures and raise after validation completes.
    EMPTY
        Replace invalid values with empty values.
    FIELDS
        Validate only the table's column structure.
    """

    # ==================================================================================
    DISABLED = "disabled"
    EAGER = "eager"
    LAZY = "lazy"
    EMPTY = "empty"
    FIELDS = "fields"
