# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import numpy as np

from pysurv.adjustment.adjustment_results import AdjustmentResults
from pysurv.basic.basic import from_rad

from .adjustment_report import AdjustmentReport


class Report(AdjustmentReport):
    """Class for representing adjustment results."""

    @property
    def general_information(self):
        content = [
            self.calculation_information,
            self.adjustment_method_information,
            self.network_specification,
            self.residual_sigma_information,
        ]
        return "\n\n".join(content)

    def __str__(self):
        if not self._results:
            return ""

        report_content = [
            self._title,
            self.general_information,
            self.controls_information,
            self.observations_information,
        ]

        return "\n\n".join(report_content)

    def to_txt(self, filename: str) -> None:
        controls_information = self._controls_information_header + "\n"
        controls_information += self._controls_information.to_string()

        observations_information = self._observations_information_header + "\n"
        observations_information += self._observations_information.to_string()

        report_content = [
            self._title,
            self.general_information,
            controls_information,
            observations_information,
        ]

        with open(filename, "w") as report_file:
            report_file.write("\n\n".join(report_content))

    def _prepare_controls_information_table(self):
        controls_table = (
            self._results.approx_coordinates.join(
                self._results.coord_corrections,
            )
            .join(
                self._results.dataset.controls.coordinates,
            )
            .join(
                self._results.adjusted_coordinate_sigmas,
            )
        )

        error_ellipses = self._results.coordinate_error_elipses.copy()
        error_ellipses["phi"] = from_rad(
            error_ellipses["phi"], self._results.dataset.measurements.angle_unit
        )
        controls_table = controls_table.join(error_ellipses).join(
            self._results.normalized_corrections
        )

        return controls_table

    def _prepare_observations_information_table(self):
        approx_observation = self._dataset.measurements_view.drop(
            columns=self._dataset.measurements.sigma_columns
        ).rename(
            columns={col: f"{col}_0" for col in self._dataset.measurements.columns},
        )

        obs_residuals = self._results.obs_residuals.copy()
        ang_columns = [
            f"v{col}" for col in self._dataset.measurements.angular_measurement_columns
        ]
        obs_residuals[ang_columns] = from_rad(
            obs_residuals[ang_columns], self._dataset.measurements.angle_unit
        )

        obs_table = approx_observation.join(obs_residuals)

        obs_values = self._results.adjusted_observation_values.copy()
        angular_measurement_columns = (
            self._dataset.measurements.angular_measurement_columns
        )
        obs_values[angular_measurement_columns] = from_rad(
            obs_values[angular_measurement_columns],
            self._dataset.measurements.angle_unit,
        )
        obs_table = obs_table.join(obs_values)

        obs_sigmas = self._results.adjusted_observation_sigmas.copy()
        angular_sigma_columns = self._dataset.measurements.angular_sigma_columns
        obs_sigmas[angular_sigma_columns] = from_rad(
            obs_sigmas[angular_sigma_columns], self._dataset.measurements.angle_unit
        )
        obs_table = obs_table.join(obs_sigmas)

        obs_table = obs_table.join(self._results.normalized_residuals)

        return obs_table

    def _format_tuning_constants(self, tuning_constants: dict | None):
        if tuning_constants is None:
            return
        return {
            key: float(np.round(value, 5)) for key, value in tuning_constants.items()
        }
