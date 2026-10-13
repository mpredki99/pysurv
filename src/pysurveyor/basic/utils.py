"""Shared utilities for numerical operations."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

# ======================================================================================
FloatArray = NDArray[np.float64]


# ======================================================================================
def is_scalar(value: object) -> bool:
    """Return whether value is a supported scalar numeric type."""
    return isinstance(value, (int, float, np.number)) or (
        isinstance(value, np.ndarray) and value.ndim == 0
    )


# ======================================================================================
def as_scalar(value: object) -> float:
    """Convert value to a scalar."""
    return float(value)


# ======================================================================================
def as_array(value: object) -> FloatArray:
    """Convert value to a NumPy array with float64 precision."""
    return np.asarray(value, dtype=np.float64)


# ======================================================================================
def normalize_numeric_type(value: object) -> float | FloatArray:
    """Normalize the numeric type of value to a supported scalar or array type."""
    if is_scalar(value):
        return as_scalar(value)
    return as_array(value)
