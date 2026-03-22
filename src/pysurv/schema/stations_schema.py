# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from pysurv.schema.utils import SchemaValidationMode, get_models_dir

from .strict_schema import StrictSchema


class StationsSchema(StrictSchema):

    def __init__(
        self,
        model: pd.DataFrame | None = None,
        validation_mode: SchemaValidationMode | str = SchemaValidationMode.LAZY,
        station_attribute_columns: Iterable[str] = ["stn_h", "stn_sh", "rz", "srz"],
    ) -> None:
        if model is None:
            model_path = self._get_stations_model_file_path()
            model = pd.read_csv(model_path)

        self.station_attribute_columns = pd.Index(station_attribute_columns)

        super().__init__(model, validation_mode=validation_mode)

    # ----------------------------------------------------------------------------------
    # Internal helpers
    # ----------------------------------------------------------------------------------
    def _get_stations_model_file_path(self) -> Path:
        models_dir = get_models_dir()
        model_path = models_dir / "stations_model.csv"

        if not model_path.is_file():
            raise FileNotFoundError(f"Stations model file not found in: {model_path}")
        return model_path
