# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def valid_data():
    return pd.DataFrame.from_dict(
        {
            "row_1": ["sample", "text field", 100, 0],
            "row_2": [0, "text value", 200, 100],
            "row_3": [True, "description", 300, None],
            "row_4": [None, 100, None, "400"],
        },
        columns=["unvalidated", "text", "distance", "angle"],
        orient="index",
    ).rename_axis("field")


@pytest.fixture
def invalid_data():
    return pd.DataFrame.from_dict(
        {
            "row_1": ["sample", "text field", 100, 0],
            "row_2": [0, None, 200, 100],
            "row_3": [True, "description", "text", None],
            "row_4": [None, 100, None, "400"],
        },
        columns=["unvalidated", "text", "distance", "angle"],
        orient="index",
    ).rename_axis("field")


@pytest.fixture
def missing_column_data():
    return pd.DataFrame.from_dict(
        {
            "row_1": ["text field", 100, 0],
            "row_2": [None, 200, 100],
            "row_3": ["description", "text", None],
            "row_4": [100, None, "400"],
        },
        columns=["text", "distance", "angle"],
        orient="index",
    ).rename_axis("field")


@pytest.fixture
def model_file_path():
    return Path(__file__).parent / "test_model.csv"
