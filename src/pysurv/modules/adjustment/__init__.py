# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

"""
adjustment submodule implements functions to define and solve a system of observation equations by reducing
them to a system of normal equations and solving this system by the least square method.
"""
from . import computations
from . import matrices
from . import observation_equations
from . import robust
