# Coding: UTF-8

# Copyright (C) 2025 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
Basic surveying functions.

This module provides fundamental surveying functions
commonly used in surveying applications.

Functions
---------
to_rad : Convert angles from degrees or gradians to radians
from_rad : Convert angles from radians to degrees or gradians
azimuth : Calculate azimuth from coordinate differences

Notes
-----
All functions support both scalar and array inputs due to implementation
in numpy.
"""

from ._constants import RHO
from .basic import azimuth, from_rad, to_rad

__all__ = ["RHO", "azimuth", "from_rad", "to_rad"]
