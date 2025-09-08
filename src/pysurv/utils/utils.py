# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from __future__ import annotations

from functools import cached_property, wraps
from typing import Any, Callable, Type

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


class refreshable_property:
    """
    Descriptor that creates a cachable property that can be refreshed
    if `refresh` condition is True. If reset all is True, it will reset
    values of all cached properties that are stored in the object
    instance (including ordinary cached properties).
    """

    def __init__(
        self,
        func: Callable | None = None,
        refresh: bool | str | Callable = True,
        reset_all: bool = False,
    ) -> None:
        self.func = func
        self.refresh = refresh
        self.reset_all = reset_all

        if func:
            wraps(func)(self)

    def __call__(self, func: Callable) -> refreshable_property:
        """Enable decorator usage without parentheses."""
        self.func = func
        wraps(func)(self)
        return self

    def __set_name__(self, owner: Type[Any], name: str) -> None:
        """Store the property name when class is created."""
        self.name = name

    def __set__(self, instance: Any, value: Any) -> None:
        """Assign a value to the property cache."""
        instance.__dict__[self.name] = value

    def __get__(self, instance: Any, owner: Type[Any] | None = None) -> Any:
        """Return the cached value or recompute it if needed."""
        if instance is None:
            return self
        return self._get_value(instance)

    def __delete__(self, instance: Any) -> None:
        """Delete the cached property value."""
        instance.__dict__.pop(self.name, None)

    def _get_value(self, instance: Any) -> Any:
        """Retrieve the property value, refreshing it if the condition is met."""
        if self._do_refresh(instance) or self.name not in instance.__dict__:

            if self.reset_all:
                reset_object_cache(instance, deep=True)

            instance.__dict__[self.name] = self.func(instance)
        return instance.__dict__[self.name]

    def _do_refresh(self, instance: Any) -> bool:
        """Determine if the property should be refreshed."""
        if isinstance(self.refresh, bool):
            return self.refresh
        if isinstance(self.refresh, str):
            return getattr(instance, self.refresh)
        if callable(self.refresh):
            return self.refresh(instance)
        return True


def reset_object_cache(obj: object, deep=True) -> None:
    """Reset object's cached properties values including cached refreshable properties."""
    mro = obj.__class__.__mro__ if deep else [obj.__class__]
    for cls in mro:
        for name, attr in cls.__dict__.items():
            if isinstance(attr, (cached_property, refreshable_property)):
                obj.__dict__.pop(name, None)


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
