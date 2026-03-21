# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from functools import cached_property
from pathlib import Path

import pandas as pd

from pysurv.schema.utils import SchemaValidationMode, get_models_dir

from .strict_schema import StrictSchema


class ControlsSchema(StrictSchema):

    def __init__(
        self,
        model: pd.DataFrame | None = None,
        validation_mode: SchemaValidationMode | str = SchemaValidationMode.LAZY,
    ) -> None:
        if model is None:
            model_path = self._get_controls_model_file_path()
            model = pd.read_csv(model_path)

        super().__init__(model, validation_mode=validation_mode)

    # ----------------------------------------------------------------------------------
    # Properties
    # ----------------------------------------------------------------------------------
    @cached_property
    def coordinate_columns(self) -> pd.Index:
        return pd.Index(["x", "y", "z"])

    @cached_property
    def coordinate_sigma_columns(self) -> pd.Index:
        return pd.Index(["sx", "sy", "sz"])

    # ----------------------------------------------------------------------------------
    # Internal helpers
    # ----------------------------------------------------------------------------------
    def _get_controls_model_file_path(self) -> Path:
        models_dir = get_models_dir()
        model_path = models_dir / "controls_model.csv"

        if not model_path.is_file():
            raise FileNotFoundError(f"Controls model file not found in: {model_path}")
        return model_path
