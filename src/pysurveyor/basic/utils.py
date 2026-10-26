# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Shared utilities for numerical operations."""

from __future__ import annotations

from numbers import Number
from typing import SupportsFloat, TypeGuard

import numpy as np
from numpy.typing import NDArray

# ======================================================================================
FloatArray = NDArray[np.float64]


# ======================================================================================
def is_scalar(value: object) -> TypeGuard[Number | np.number | np.ndarray]:
    """Return whether value is a supported scalar numeric type."""
    return isinstance(value, (Number, np.number)) or (
        isinstance(value, np.ndarray) and value.ndim == 0
    )


# ======================================================================================
def as_scalar(value: Number | np.number | np.ndarray) -> float:
    """Convert value to a scalar."""
    if isinstance(value, np.ndarray):
        if value.ndim != 0:
            raise ValueError("Expected a zero-dimensional array.")
        return float(value.item())
    # ----------------------------------------------------------------------------------
    if isinstance(value, np.number):
        return float(value)
    # ----------------------------------------------------------------------------------
    if isinstance(value, SupportsFloat):
        return float(value)
    # ----------------------------------------------------------------------------------
    raise TypeError(f"Cannot convert {type(value).__name__} to float.")


# ======================================================================================
def as_array(value: object) -> FloatArray:
    """Convert value to a NumPy array with float64 precision."""
    return np.asarray(value, dtype=np.float64)


# ======================================================================================
def normalize_numeric_type(value: object) -> float | FloatArray:
    """Normalize the numeric type of value to a supported scalar or array type."""
    if is_scalar(value):
        return as_scalar(value)
    # ----------------------------------------------------------------------------------
    return as_array(value)
