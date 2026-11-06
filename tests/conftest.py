# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def tmp_text_file(tmp_path):
    """Create a temporary text file."""

    def _create(
        filename: str = "test.txt", data: pd.DataFrame | str | None = None
    ) -> Path:
        filepath = tmp_path / filename

        # ------------------------------------------------------------------------------
        if data is None:
            filepath.write_text("")

        # ------------------------------------------------------------------------------
        if isinstance(data, pd.DataFrame):
            data.to_csv(filepath, index=False)

        # ------------------------------------------------------------------------------
        if isinstance(data, str):
            filepath.write_text(data)

        # ------------------------------------------------------------------------------
        return filepath

    # ----------------------------------------------------------------------------------
    return _create
