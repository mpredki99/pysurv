# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from enum import StrEnum
from pathlib import Path


class SchemaValidationMode(StrEnum):
    """Defines behavior of `validate` method in `PySurvSchema`."""

    DISABLED = "disabled"  # Skip validation
    EAGER = "eager"  # Raise ValidationError on first failure
    LAZY = "lazy"  # Collect errors, raise after validation
    COERCE = "coerce"  # Try to coerce values into correct type
    FIELDS = "fields"  # Validate only columns structure


def get_models_dir() -> Path:
    """Return path object containing the `models` directory."""
    models_dir = Path(__file__).parent.parent / "models"

    if not models_dir.is_dir():
        raise FileNotFoundError(f"Models directory does not exist: {models_dir}")
    return models_dir
