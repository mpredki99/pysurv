# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from functools import cached_property
from itertools import chain

import pandas as pd

from pysurv.schema.utils import SchemaValidationMode, get_models_dir

from .strict_schema import StrictSchema


class MeasurementsSchema(StrictSchema):

    def __init__(
        self,
        model: pd.DataFrame | None = None,
        validation_mode: SchemaValidationMode | str = SchemaValidationMode.LAZY,
    ) -> None:
        if model is None:
            model_path = self._get_measurements_model_file_path()
            model = pd.read_csv(model_path)

        super().__init__(model, validation_mode=validation_mode)

    # ----------------------------------------------------------------------------------
    # Properties
    # ----------------------------------------------------------------------------------
    @cached_property
    def target_columns(self) -> pd.Index:
        return pd.Index(["trg_id", "trg_h", "trg_sh", "trg_ctr", "trg_cst", "trg_scst"])

    @cached_property
    def linear_measurements_columns(self):
        return pd.Index(["sd", "hd", "vd", "dx", "dy", "dz"])

    @cached_property
    def linear_measurements_sigma_columns(self):
        return pd.Index(["ssd", "shd", "svd", "sdx", "sdy", "sdz"])

    @cached_property
    def angular_measurements_columns(self):
        return pd.Index(["a", "hz", "vz", "vh"])

    @cached_property
    def angular_measurements_sigma_columns(self):
        return pd.Index(["sa", "shz", "svz", "svh"])

    @cached_property
    def measurements_columns(self):
        return self.linear_measurements_columns.union(self.angular_measurements_columns)

    @cached_property
    def measurements_sigma_columns(self):
        return self.linear_measurements_sigma_columns.union(
            self.angular_measurements_sigma_columns
        )

    @cached_property
    def linear_columns(self):
        return pd.Index(
            chain.from_iterable(
                zip(
                    self.linear_measurements_columns,
                    self.linear_measurements_sigma_columns,
                )
            )
        )

    @cached_property
    def angular_columns(self):
        return pd.Index(
            chain.from_iterable(
                zip(
                    self.angular_measurements_columns,
                    self.angular_measurements_sigma_columns,
                )
            )
        )

    # ----------------------------------------------------------------------------------
    # Internal helpers
    # ----------------------------------------------------------------------------------
    def _get_measurements_model_file_path(self):
        models_dir = get_models_dir()
        model_path = models_dir / "measurements_model.csv"

        if not model_path.is_file():
            raise FileNotFoundError(
                f"Measurements model file does not exist: {model_path}"
            )
        return model_path
