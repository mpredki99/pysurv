# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import re
from abc import ABC, abstractmethod
from collections.abc import Iterable
from typing import Any, Callable

import pandas as pd

from pysurv.exceptions import ValidationError

from .constants import COMMENT
from .utils import PySurvValidatorMode


class PySurvValidator(ABC):
    def __init__(
        self,
        *,
        ignore: Iterable[Any] | Any = (pd.NA, COMMENT),
        mode: PySurvValidatorMode | str = PySurvValidatorMode.RAISE,
    ) -> None:
        """Configure the validator."""
        super().__init__()
        self.ignore = ignore
        self.mode = mode

    # ----------------------------------------------------------------------------------
    # Dunder methods
    # ----------------------------------------------------------------------------------
    @abstractmethod
    def __call__(
        self,
        data: pd.Series,
    ) -> tuple[pd.Series, pd.Series]:
        """Handle validation logic."""
        ...

    def __and__(self, other: "PySurvValidator") -> "AndValidator":
        """Allow chaining validators."""
        if not isinstance(other, PySurvValidator):
            raise TypeError("Cannot combine PySurvValidator with non-PySurvValidator")
        return AndValidator(self, other)

    def __str__(self) -> str:
        return self.__class__.__name__

    def __repr__(self) -> str:
        """String representation with parameter names (leading '_' are stripped)."""
        attrs = ", ".join(
            f"{key.lstrip('_')}={value!r}" for key, value in self.__dict__.items()
        )
        return f"{self.__class__.__name__}({attrs})"

    def __contains__(self, value: "PySurvValidator") -> bool:
        name = getattr(value, "__name__", None)
        if name:
            return self.__class__.__name__ == name
        return isinstance(self, type(value))

    # ----------------------------------------------------------------------------------
    # Properties
    # ----------------------------------------------------------------------------------
    @property
    def ignore(self) -> dict[str, Any]:
        """Return current ignore settings as a dict."""
        ignored = {"na": self._ignore_na}

        if self._ignore_literals:
            ignored["literals"] = set(self._ignore_literals)

        if self._ignore_regexes:
            ignored["regexes"] = tuple(self._ignore_regexes)

        if self._ignore_callables:
            ignored["callables"] = tuple(self._ignore_callables)

        return ignored

    @ignore.setter
    def ignore(self, values: Iterable[Any] | Any) -> None:
        """
        Configure which values should be ignored during validation.

        This method accepts a single value or an iterable of values to specify which entries
        in the data should be excluded from validation checks:

        - An empty iterable means no values will be ignored
        - `pd.NA` and `None` are treated as missing values to be ignored
        - Strings are treated as scalar literal values to ignore
        - Regular expressions must be compiled with `re.compile`;
          matching these will be ignored
        - Callables have to return boolean values;
          will be applied to non-empty values;
          if they return True, the value is ignored
        """
        if isinstance(values, Iterable) and not isinstance(values, (str, bytes)):
            raw_values = values
        else:
            raw_values = [values]

        ignore_na = False
        literals: set[Any] = set()
        regexes: list[re.Pattern] = []
        callables: list[Callable[[Any], bool]] = []

        for v in raw_values:
            if v is None or v is pd.NA:
                ignore_na = True

            elif isinstance(v, re.Pattern):
                regexes.append(v)

            elif callable(v):
                callables.append(v)

            else:
                literals.add(v)

        self._ignore_na = ignore_na
        self._ignore_literals = literals
        self._ignore_regexes = regexes
        self._ignore_callables = callables

    @property
    def mode(self) -> PySurvValidatorMode:
        return self._mode

    @mode.setter
    def mode(self, value: PySurvValidatorMode | str) -> None:
        """Set PySurvValidatorMode or create it from string."""
        if isinstance(value, PySurvValidatorMode):
            self._mode = value
        else:
            self._mode = PySurvValidatorMode(value)

    # ----------------------------------------------------------------------------------
    # Internal helpers
    # ----------------------------------------------------------------------------------
    def _ignored_mask(self, data: pd.Series) -> pd.Series:
        """Determine ignored values mask."""
        mask = pd.Series(False, index=data.index)

        # Missing values
        if self._ignore_na:
            mask |= data.isna()

        # Literal values
        if self._ignore_literals:
            mask |= data.isin(self._ignore_literals)

        # Regex patterns
        if self._ignore_regexes:
            str_values = data.astype("string")

            for pattern in self._ignore_regexes:
                mask |= str_values.str.match(pattern, na=False)

        # Callables
        if self._ignore_callables:
            for func in self._ignore_callables:
                mask |= data.apply(lambda x: False if pd.isna(x) else func(x))

        return mask

    def _return_validation_result(
        self,
        mask: pd.Series,
        data: pd.Series,
        error_message: str = "Invalid values found:",
    ) -> tuple[pd.Series, pd.Series]:
        """Return validation results or raise validation error."""
        if self.mode == PySurvValidatorMode.RAISE and not mask.all():
            self._raise_validation_error(data[~mask], error_message)
        else:
            return mask, data

    @staticmethod
    def _raise_validation_error(invalid_values: pd.Series, message: str) -> None:
        """Raise formatted error message."""
        invalid_values = invalid_values.rename("value").rename_axis("row").reset_index()
        raise ValidationError(f"{message}\n{invalid_values.to_string(index=False)}")

    @staticmethod
    def _return_empty(
        dtype: type | str, index: pd.Index
    ) -> tuple[pd.Series, pd.Series]:
        """Return empty pandas Series with proper dtype and index."""
        return pd.Series(dtype=bool, index=index), pd.Series(dtype=dtype, index=index)


# --------------------------------------------------------------------------------------
#          AndValidator
# --------------------------------------------------------------------------------------
class AndValidator(PySurvValidator):
    def __init__(
        self,
        left: PySurvValidator,
        right: PySurvValidator,
        mode: PySurvValidatorMode = PySurvValidatorMode.RAISE,
    ) -> None:
        self.left = left
        self.right = right
        self.mode = mode

    # ----------------------------------------------------------------------------------
    # Dunder methods
    # ----------------------------------------------------------------------------------
    def __call__(self, data: pd.Series) -> tuple[pd.Series, pd.Series]:
        """Validate data using both validators."""
        mask_1, validated = self.left(data)
        mask_2, validated = self.right(validated)

        valid_mask = mask_1 & mask_2

        return self._return_validation_result(valid_mask, validated)

    def __str__(self) -> str:
        return f"{self.left} & {self.right}"

    def __repr__(self) -> str:
        return f"{self.left!r} & {self.right!r}"

    def __contains__(self, value: PySurvValidator) -> bool:
        return value in self.left or value in self.right

    def __iter__(self) -> Iterable:
        for validator in (self.left, self.right):
            if isinstance(validator, AndValidator):
                yield from validator
            else:
                yield validator

    def __len__(self) -> int:
        left = len(self.left) if isinstance(self.left, AndValidator) else 1
        right = len(self.right) if isinstance(self.right, AndValidator) else 1
        return left + right

    def __getitem__(self, key):
        validators = tuple(self)

        if isinstance(key, int):
            return validators[key]

        if isinstance(key, slice):
            subset = validators[key]
            result = subset[0]

            for validator in subset[1:]:
                result = result & validator

            return result

    # ----------------------------------------------------------------------------------
    # Properties
    # ----------------------------------------------------------------------------------
    @property
    def ignore(self):
        raise AttributeError(
            "AndValidator does not support `ignore`. "
            "Configure ignore on individual left/right validators."
        )

    @ignore.setter
    def ignore(self, value):
        raise AttributeError(
            "Cannot set `ignore` on AndValidator. "
            "Set ignore on left/right validators instead."
        )

    @property
    def left(self) -> PySurvValidator:
        return self._left

    @left.setter
    def left(self, validator: PySurvValidator) -> None:
        validator.mode = PySurvValidatorMode.RETURN
        self._left = validator

    @property
    def right(self) -> PySurvValidator:
        return self._right

    @right.setter
    def right(self, validator: PySurvValidator) -> None:
        validator.mode = PySurvValidatorMode.RETURN
        self._right = validator
