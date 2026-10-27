# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.backend.backend import DataFrame
from pysurveyor.schema.model.control_points_model import get_control_points_model_data
from pysurveyor.schema.model.instruments_model import get_instruments_model_data
from pysurveyor.schema.model.measurements_model import get_measurements_model_data
from pysurveyor.schema.model.points_model import get_points_model_data
from pysurveyor.schema.model.stations_model import get_stations_model_data
from pysurveyor.schema.model.targets_model import get_targets_model_data


# ======================================================================================
# Schema models factory class
# ======================================================================================
class Models:
    """Class representing the schema models for PySurveyor."""

    # ==================================================================================
    # Properties for lazy loading of models
    # ==================================================================================
    @property
    def points(self) -> DataFrame:
        data = get_points_model_data()
        return self._build_model(data)

    # ==================================================================================
    @property
    def control_points(self) -> DataFrame:
        data = get_control_points_model_data()
        return self._build_model(data)

    # ==================================================================================
    @property
    def stations(self) -> DataFrame:
        data = get_stations_model_data()
        return self._build_model(data)

    # ==================================================================================
    @property
    def instruments(self) -> DataFrame:
        data = get_instruments_model_data()
        return self._build_model(data)

    # ==================================================================================
    @property
    def targets(self) -> DataFrame:
        data = get_targets_model_data()
        return self._build_model(data)

    # ==================================================================================
    @property
    def measurements(self) -> DataFrame:
        data = get_measurements_model_data()
        return self._build_model(data)

    # ==================================================================================
    # Private method to build a DataFrame from model data
    # ==================================================================================
    @staticmethod
    def _build_model(model_data) -> DataFrame:
        """Builds a DataFrame from the provided model data dictionary."""

        model = DataFrame.from_dict(
            model_data,
            columns=["required", "unit", "validator", "description"],
            orient="index",
        )
        # ------------------------------------------------------------------------------
        model.index.name = "field"

        # ------------------------------------------------------------------------------
        return model


# ======================================================================================
# Entry point for schema models
# ======================================================================================
models = Models()
