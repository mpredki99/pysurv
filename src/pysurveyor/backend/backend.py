"""Pandas-backed types and functions used throughout PySurveyor.

This module centralizes the small set of pandas types and functions that are
part of the PySurveyor typing/backend interface.
"""

from __future__ import annotations

import pandas as pd


DataFrame = pd.DataFrame

Index = pd.Index

read_csv = pd.read_csv
