# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from dataclasses import asdict, dataclass
from typing import Any, Callable


@dataclass
class CSVConfig:
    delimiter: str | None = None
    decimal: str = "."
    skiprows: int | Callable | None = None
    nrows: int | None = None

    @property
    def kwargs(self) -> dict[str, Any]:
        return {k: v for k, v in asdict(self).items() if v is not None}
