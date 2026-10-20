# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Pandas-backed types and functions used throughout PySurveyor."""

from __future__ import annotations

import pandas as pd

DataFrame = pd.DataFrame
Series = pd.Series
Index = pd.Index
NA = pd.NA

isna = pd.isna
to_numeric = pd.to_numeric

read_csv = pd.read_csv
