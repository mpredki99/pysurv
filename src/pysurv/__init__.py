# Coding: UTF-8

"""
************************************************************************************************************************

    Ordinary, weighted, and robust least squares adjustment of surveying control networks.
    Copyright (C) 2024, Michal Predki

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

************************************************************************************************************************

pysurv is a Python package for adjusting surveying control networks.
The package supports importing data from CSV files and performing
ordinary, weighted, and robust least squares adjustment.
It also allows the free adjustment approach to be combined with
ordinary, weighted, and robust methods.
Additionally, you can mix these methods when adjusting
observations and reference points.
After completing the calculations, a detailed HTML report with
the adjustment results can be generated.
"""
__author__ = 'Michał Prędki'
__date__ = 'September 2024'

from .project_menager import Project
from . import customizations
from . import modules
