# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

"""Pandas-backed types and functions used throughout PySurveyor.

This module centralizes the small set of pandas types and functions that are
part of the PySurveyor typing/backend interface.
"""

from __future__ import annotations

import pandas as pd

DataFrame = pd.DataFrame

Index = pd.Index

read_csv = pd.read_csv
