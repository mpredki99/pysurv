# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from functools import cached_property, wraps
from typing import Callable

import numpy as np


def inf_to_zero(func):
    """Decorator that turns infinite values to 0 (robust methods have limit 0 with v -> inf)."""

    @wraps(func)
    def wrapper(v, *args, **kwargs):
        is_finite_mask = np.isfinite(v)
        coeff = np.nan_to_num(v, nan=np.nan, neginf=0, posinf=0)
        coeff[is_finite_mask] = func(v[is_finite_mask], *args, **kwargs)
        return coeff

    return wrapper


def apply_where(
    arg: np.ndarray,
    cond: np.ndarray,
    func_true: Callable[[np.ndarray], np.ndarray],
    func_false: Callable[[np.ndarray], np.ndarray],
) -> np.ndarray:
    """
    Apply `func_true` to elements of `arg` where `cond` is True,
    and `func_false` to elements where `cond` is False.

    Returns a new array with the same shape as `arg`.
    """
    result = np.empty(arg.shape)

    if np.any(cond):
        result[cond] = func_true(arg[cond])
    if np.any(~cond):
        result[~cond] = func_false(arg[~cond])

    return result


def reset_object_cache(object: object, deep=True) -> None:
    """Reset object's cached properties values."""
    mro = object.__class__.__mro__ if deep else [object.__class__]
    for cls in mro:
        for name, attr in cls.__dict__.items():
            if isinstance(attr, cached_property):
                object.__dict__.pop(name, None)
