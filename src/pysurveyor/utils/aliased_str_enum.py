# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from enum import StrEnum
from typing import Self


class AliasedStrEnum(StrEnum):
    """Class for creating string enumerations with aliases.

    This class extends the built-in StrEnum class to allow for the definition
    of string enumerations with multiple aliases for each member. It provides
    a mechanism to look up enumeration members using either their canonical
    value or any of their defined aliases. It also supports case-insensitive
    lookups if the __to_lower__ attribute is set to True.
    """

    __to_lower__ = False
    _aliases_: tuple[str, ...]

    def __new__(cls, value: str, *aliases: str):
        self = str.__new__(cls, value)
        self._value_ = value
        self._aliases_ = aliases

        return self

    @classmethod
    def _missing_(cls, value) -> Self | None:
        if not isinstance(value, str):
            return None

        normalized_value = value.lower() if cls.__to_lower__ else value

        for member in cls:
            if member.value == normalized_value:
                return member

            if any(alias == normalized_value for alias in member._aliases_):
                return member

        return None
