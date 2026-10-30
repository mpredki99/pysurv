# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from abc import abstractmethod
from collections.abc import Iterable, Iterator
import re
from typing import Any, Callable

from pysurveyor.abc.pysurveyor_abc import PySurveyorABC
from pysurveyor.backend.backend import NA, Index, Series, isna
from pysurveyor.exceptions.exceptions import ValidationError
from pysurveyor.validators.utils.constants import COMMENT
from pysurveyor.validators.utils.validator_mode import ValidatorMode


# ======================================================================================
# Abstract base class for PySurveyor validators
# ======================================================================================
class Validator(PySurveyorABC):
    """Base class for validators that check a Series and report row-level results.

    Subclasses implement ``__call__`` and return a boolean validity mask together
    with the data. Validators can be combined with ``&``.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    # Dunder methods
    # ==================================================================================
    def __init__(
        self,
        *,
        ignore: Iterable[Any] | Any = (NA, COMMENT),
        mode: ValidatorMode | str = ValidatorMode.RAISE,
    ) -> None:
        """Initialize ignore rules and choose how validation results are handled."""
        super().__init__()
        # ------------------------------------------------------------------------------
        self._ignore = ignore
        self._ignore_na = False
        self._ignore_literals: set[Any] = set()
        self._ignore_regexes: list[re.Pattern] = []
        self._ignore_callables: list[Callable[[Any], bool]] = []
        # ------------------------------------------------------------------------------
        self.ignore = ignore
        self.mode = mode

    # ==================================================================================
    @abstractmethod
    def __call__(
        self,
        data: Series,
    ) -> tuple[Series, Series]:
        """Validate ``data`` and return ``(validity_mask, data)``.

        The mask is indexed like ``data`` and contain ``True`` for rows
        considered valid. Implementations may return transformed data.
        """
        pass

    # ==================================================================================
    def __and__(self, other: Validator) -> AndValidator:
        """Allow chaining validators using logical AND."""
        if not isinstance(other, Validator):
            raise TypeError("Cannot combine Validator with non-Validator")
        return AndValidator(self, other)

    # ==================================================================================
    def __str__(self) -> str:
        """Return the validator class name."""
        return self.__class__.__name__

    # ==================================================================================
    def __repr__(self) -> str:
        """Return a constructor-like representation of the validator's attributes."""
        attrs = ", ".join(
            f"{key.lstrip('_')}={value!r}"
            for key, value in self.__dict__.items()
            if key
            not in (
                "_ignore_na",
                "_ignore_literals",
                "_ignore_regexes",
                "_ignore_callables",
            )
        )
        return f"{self.__class__.__name__}({attrs})"

    # ==================================================================================
    def __contains__(self, value: Validator) -> bool:
        """Check whether ``value`` identifies this validator by name or type."""
        # ------------------------------------------------------------------------------
        name = getattr(value, "__name__", None)
        if name is not None:
            return self.__class__.__name__ == name
        # ------------------------------------------------------------------------------
        return isinstance(self, type(value))

    # ==================================================================================
    # Properties
    # ==================================================================================
    @property
    def ignore(self) -> Iterable[Any] | Any:
        """Return the configured ignore rules."""
        return self._ignore

    @property
    def ignore_dict(self) -> dict[str, Any]:
        """Expose configured ignore rules grouped by rule type."""
        ignored = {}

        # ------------------------------------------------------------------------------
        ignored["empty"] = self._ignore_na

        # ------------------------------------------------------------------------------
        if self._ignore_literals:
            ignored["literals"] = set(self._ignore_literals)

        # ------------------------------------------------------------------------------
        if self._ignore_regexes:
            ignored["regexes"] = list(self._ignore_regexes)

        # ------------------------------------------------------------------------------
        if self._ignore_callables:
            ignored["callables"] = list(self._ignore_callables)

        # ------------------------------------------------------------------------------
        return ignored

    @ignore.setter
    def ignore(self, values: Iterable[Any] | Any) -> None:
        """Configure values that should be skipped by validation.

        ``values`` may be a scalar or an iterable of values to ignore.
        The following types are supported:

        - An empty iterable means no values will be ignored
        - if ``NA`` or ``None`` is included, missing values will be ignored
        - Strings are treated as scalar literal values to ignore
        - Regular expressions compiled with `re.compile`
        - Callables that return boolean values - if True, then the value is ignored
        """
        raw_values = (
            values
            if isinstance(values, Iterable) and not isinstance(values, (str, bytes))
            else (values,)
        )

        # ------------------------------------------------------------------------------
        ignore_na = False
        literals = set()
        regexes = []
        callables = []

        # ------------------------------------------------------------------------------
        for value in raw_values:
            if value is None or value is NA:
                ignore_na = True
            # --------------------------------------------------------------------------
            elif isinstance(value, re.Pattern):
                regexes.append(value)
            # --------------------------------------------------------------------------
            elif callable(value):
                callables.append(value)
            # --------------------------------------------------------------------------
            else:
                literals.add(value)

        # ------------------------------------------------------------------------------
        self._ignore = values
        self._ignore_na = ignore_na
        self._ignore_literals = literals
        self._ignore_regexes = regexes
        self._ignore_callables = callables

    # ==================================================================================
    @property
    def mode(self) -> ValidatorMode:
        """Return the current validation mode."""
        return self._mode

    @mode.setter
    def mode(self, value: ValidatorMode | str) -> None:
        """Set the validation mode, converting strings to ``ValidatorMode``."""
        self._mode = value if isinstance(value, ValidatorMode) else ValidatorMode(value)

    # ==================================================================================
    # Private methods
    # ==================================================================================
    def _ignored_mask(self, data: Series) -> Series:
        """Build a boolean mask whose true entries should be excluded from checks."""
        mask = Series(False, index=data.index)

        # Missing values ---------------------------------------------------------------
        if self._ignore_na:
            mask |= data.isna()

        # Literal values ---------------------------------------------------------------
        if self._ignore_literals:
            mask |= data.isin(self._ignore_literals)

        # Regex patterns ---------------------------------------------------------------
        if self._ignore_regexes:
            string_values = data.astype("string")
            for pattern in self._ignore_regexes:
                mask |= string_values.str.match(pattern, na=False)

        # Callables --------------------------------------------------------------------
        if self._ignore_callables:
            for predicate in self._ignore_callables:
                mask |= data.apply(
                    lambda value, predicate=predicate: (
                        False if isna(value) else predicate(value)
                    )
                )

        # ------------------------------------------------------------------------------
        return mask

    # ==================================================================================
    def _return_validation_result(
        self,
        mask: Series,
        data: Series,
        error_message: str = "Invalid values found:",
    ) -> tuple[Series, Series]:
        """Return the mask and data, or raise if RAISE mode finds invalid rows.

        Parameters
        ----------
        mask:
            Boolean Series with ``True`` for valid rows.
        data:
            Data associated with the validation result.
        error_message:
            Prefix used when formatting a validation exception.
        """
        if self.mode == ValidatorMode.RAISE and not mask.all():
            self._raise_validation_error(data[~mask], error_message)
        return mask, data

    # ==================================================================================
    @staticmethod
    def _raise_validation_error(invalid_values: Series, message: str) -> None:
        """Raise ``ValidationError`` with invalid values and their row labels."""
        formatted = invalid_values.rename("value").rename_axis("row").reset_index()
        raise ValidationError(f"{message}\n{formatted.to_string(index=False)}")

    # ==================================================================================
    @staticmethod
    def _return_empty(dtype: type | str, index: Index) -> tuple[Series, Series]:
        """Create an empty boolean result Series and an empty data Series.

        Both Series use the supplied index; the data Series uses ``dtype``.
        """
        return Series(dtype=bool, index=index), Series(dtype=dtype, index=index)


# ======================================================================================
# AndValidator for combining two validators with logical AND
# ======================================================================================
class AndValidator(Validator):

    # ==================================================================================
    # Dunder methods
    # ==================================================================================
    def __init__(
        self,
        left: Validator,
        right: Validator,
        mode: ValidatorMode = ValidatorMode.RAISE,
    ) -> None:
        """Create AND validator and configure its component to return masks."""
        self.left = left
        self.right = right
        self.mode = mode

    # ==================================================================================
    def __call__(self, data: Series) -> tuple[Series, Series]:
        """Run both validators in sequence and combine their validity masks.

        The right validator receives the data returned by the left validator.
        """
        left_mask, validated = self.left(data)
        right_mask, validated = self.right(validated)
        return self._return_validation_result(left_mask & right_mask, validated)

    # ==================================================================================
    def __str__(self) -> str:
        return f"{self.left} & {self.right}"

    # ==================================================================================
    def __repr__(self) -> str:
        return f"{self.left!r} & {self.right!r}"

    # ==================================================================================
    def __contains__(self, value: Validator) -> bool:
        """Return whether either component contains the specified validator."""
        return value in self.left or value in self.right

    # ==================================================================================
    def __iter__(self) -> Iterator[Validator]:
        """Iterate over component validators, flattening nested conjunctions."""
        for validator in (self.left, self.right):
            if isinstance(validator, AndValidator):
                yield from validator
            else:
                yield validator

    # ==================================================================================
    def __len__(self) -> int:
        """Return the number of leaf validators in this conjunction."""
        return sum(
            1 if not isinstance(v, AndValidator) else len(v)
            for v in (self.left, self.right)
        )

    # ==================================================================================
    def __getitem__(self, key: int | slice) -> Validator:
        """Get a leaf validator by index or combine a slice of leaf validators.

        Integer indices follow tuple indexing semantics. A non-empty slice returns
        a conjunction of the selected validators; an empty slice raises ``IndexError``.
        """
        validators = tuple(self)

        # ------------------------------------------------------------------------------
        if isinstance(key, int):
            return validators[key]

        # ------------------------------------------------------------------------------
        if isinstance(key, slice):
            subset = validators[key]
            if not subset:
                raise IndexError("Cannot build a validator from an empty slice")

            result = subset[0]
            for validator in subset[1:]:
                result = result & validator
            return result

        # ------------------------------------------------------------------------------
        raise TypeError("Validator indices must be integers or slices")

    # ==================================================================================
    # Properties
    # ==================================================================================
    @property
    def ignore(self) -> dict[str, Any]:
        """Disallow aggregate ignore rules - configure each component independently."""
        raise AttributeError(
            "AndValidator does not support `ignore`. "
            "Configure ignore on individual left/right validators."
        )

    @ignore.setter
    def ignore(self, values: Iterable[Any] | Any) -> None:
        """Prevent ignore settings from being assigned to the aggregate validator."""
        raise AttributeError(
            "Cannot set `ignore` on AndValidator. "
            "Set ignore on left/right validators instead."
        )

    # ==================================================================================
    @property
    def left(self) -> Validator:
        """Return the left-hand component validator."""
        return self._left

    @left.setter
    def left(self, validator: Validator) -> None:
        """Set the left-hand component validator."""
        validator.mode = ValidatorMode.RETURN
        self._left = validator

    # ==================================================================================
    @property
    def right(self) -> Validator:
        """Return the right-hand component validator."""
        return self._right

    @right.setter
    def right(self, validator: Validator) -> None:
        """Set the right-hand component validator."""
        validator.mode = ValidatorMode.RETURN
        self._right = validator
