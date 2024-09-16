# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import matplotlib.pyplot as plt
import pandas as pd


def to_string(plot: plt.figure, general_info: pd.DataFrame, controls_table: pd.DataFrame,
              measurements_table: pd.DataFrame, footer: list, show_plot: bool = True) -> str:
    """
    Returns adjustment results report as the string ready to display on the screen.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - plot: (figure): Plot of sigma zero values.
    - general_info: (pd.DataFrame): DataFrame with general information about adjustment.
    - controls_table: (pd.DataFrame): DataFrame containing information about control points.
    - measurements_table: (pd.DataFrame): DataFrame containing information about measurements.
    - footer: (list): Footer of the report.
    - show_plot: (bool, optional): If True, displays the plot of sigma values. Defaults to True.

    Returns:
        str: Report of adjustment results.
    """
    if show_plot:
        plot.show()
    return f"""
LEAST SQUARES ADJUSTMENT RESULTS

GENERAL INFORMATION

{general_info.to_string(header=False)
        .replace('obs_sigma_zero', 'σ₀')
        .replace('pt_sigma_zero', 'σ₀')}


CONTROL POINTS

{controls_table.reset_index().to_string(justify='center', index=False)
    .replace('_0_', '₀')
    .replace('_s_', 'σ')
    .replace('NaN', '')}


MEASUREMENTS

{measurements_table.reset_index().to_string(justify='center', index=False)
    .replace('_0_', '₀')
    .replace('_s_', 'σ')
    .replace('NaN', '')}

{(''
  '').join(footer)}
"""
