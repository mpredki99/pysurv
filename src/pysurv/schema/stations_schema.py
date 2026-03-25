# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from collections.abc import Iterable

import pandas as pd

from pysurv.schema.utils import SchemaValidationMode

from .pysurv_table_schema import PySurvTableSchema

STATIONS_MODEL_FILE_NAME = "stations_model.csv"


class StationsSchema(PySurvTableSchema):

    def __init__(
        self,
        model: pd.DataFrame | None = None,
        validation_mode: SchemaValidationMode | str = SchemaValidationMode.LAZY,
        station_attribute_columns: Iterable[str] = ["stn_h", "stn_sh", "rz", "srz"],
    ) -> None:
        if model is None:
            model_path = self._get_model_file_path(STATIONS_MODEL_FILE_NAME)
            model = pd.read_csv(model_path)

        self.station_attribute_columns = pd.Index(station_attribute_columns)

        super().__init__(model, validation_mode=validation_mode)
