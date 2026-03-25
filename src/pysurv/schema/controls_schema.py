# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from collections.abc import Iterable

import pandas as pd

from pysurv.schema.utils import SchemaValidationMode

from .pysurv_table_schema import PySurvTableSchema

CONTROLS_MODEL_FILE_NAME = "controls_model.csv"


class ControlsSchema(PySurvTableSchema):

    def __init__(
        self,
        model: pd.DataFrame | None = None,
        validation_mode: SchemaValidationMode | str = SchemaValidationMode.LAZY,
        coordinate_columns: Iterable[str] = ["x", "y", "z"],
        coordinate_sigma_columns: Iterable[str] = ["sx", "sy", "sz"],
    ) -> None:
        if model is None:
            model_path = self._get_model_file_path(CONTROLS_MODEL_FILE_NAME)
            model = pd.read_csv(model_path)

        self.coordinate_columns = pd.Index(coordinate_columns)
        self.coordinate_sigma_columns = pd.Index(coordinate_sigma_columns)

        super().__init__(model, validation_mode=validation_mode)
