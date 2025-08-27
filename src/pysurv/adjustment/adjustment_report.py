# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from abc import ABC, abstractmethod

import numpy as np
import pandas as pd

from .adjustment_results import AdjustmentResults


class AdjustmentReport(ABC):
    def __init__(self, adjustment_results: AdjustmentResults) -> None:
        self._dataset = adjustment_results.dataset
        self._results = adjustment_results

        self._title = "----- ADJUSTMENT REPORT -----"

        self._calculation_information_header = "CALCULATION INFORMATION"
        self._calculation_information = {
            "- adjustment status:": self._results.calculation_status,
            "- number of iterations:": self._results.n_iter,
        }

        self._adjustment_method_information_header = "ADJUSTMNET METHOD INFORMATION"
        self._adjustment_method_information = {
            "- method of adjustment observations:": self._results.methods.obs_adj,
            "- observation adjustment tuning constants:": self._format_tuning_constants(
                self._results.methods.obs_tuning_constants
            ),
            "- method of free adjustment:": self._results.methods.free_adjustment,
            "- free adjustment tuning constants:": self._format_tuning_constants(
                self._results.methods.free_adj_tuning_constants
            ),
            "- free adjustment inner constraints applied:": self._results.inner_constraints,
        }

        self._network_specification_header = "NETWORK SPECIFICATION"
        self._network_specification = {
            "- number of observations:": self._results.n_measurements,
            "- number of coordinate corrections:": self._results.n_coord_corrections,
            "- number of movable tie points:": self._results.n_movable_tie_points,
            "- number of fixed tie points:": self._results.n_fixed_tie_points,
            "- number of all unknowns:": self._results.n_unknowns,
            "- number degrees of freedom:": self._results.degrees_of_freedom,
        }

        self._residual_sigma_information_header = "RESIDUAL SIGMA INFORMATION"
        self._residual_sigma_information = {
            "- residual sigma:": self._results.residual_sigma.round(5),
            "- coord correction sigma:": self._results.coord_correction_sigma.round(5),
        }

        self._controls_information_header = "CONTROLS INFORMATION"
        self._controls_information = self._prepare_controls_information_table()

        self._observations_information_header = "OBSERVATIONS INFORMATION"
        self._observations_information = self._prepare_observations_information_table()

        self._footer = [
            "Calculations made with pysurv, by Michal Predki. \n",
            "Github repository: https://github.com/mpredki99/pysurv.\n",
            "You can contact me on LinkedIn for feedback or if you discover a bug.\n",
            " ",
        ]

    def _report_part_to_string(self, part: dict) -> str:
        return "\n".join(
            f"{key} {value}" for key, value in part.items() if value is not None
        )

    @property
    def calculation_information(self):
        text = self._calculation_information_header + "\n"
        text += self._report_part_to_string(self._calculation_information)
        return text

    @property
    def adjustment_method_information(self):
        text = self._adjustment_method_information_header + "\n"
        text += self._report_part_to_string(self._adjustment_method_information)
        return text

    @property
    def network_specification(self):
        text = self._network_specification_header + "\n"
        text += self._report_part_to_string(self._network_specification)
        return text

    @property
    def residual_sigma_information(self):
        text = self._residual_sigma_information_header + "\n"
        text += self._report_part_to_string(self._residual_sigma_information)
        return text

    @property
    def controls_information(self):
        text = self._controls_information_header + "\n"
        text += str(self._controls_information.round(3))
        return text

    @property
    def observations_information(self):
        text = self._observations_information_header + "\n"
        text += str(self._observations_information.round(3))
        return text

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def to_txt(self, filename: str) -> None:
        pass

    @abstractmethod
    def _prepare_controls_information_table(self):
        pass

    @abstractmethod
    def _prepare_observations_information_table(self):
        pass

    @abstractmethod
    def _format_tuning_constants(self, tuning_constants: dict | None):
        pass
