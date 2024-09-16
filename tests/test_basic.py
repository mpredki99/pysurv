# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import pytest
import numpy as np
from pysurv.modules.basic import azimuth, to_rad, from_rad

def test_azimuth_first_quarter():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': 100, 'y': 100}
    result = azimuth(first_point, second_point)
    assert pytest.approx(result) == np.pi / 4

def test_azimuth_second_quarter():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': -100, 'y': 100}
    result = azimuth(first_point, second_point)
    assert pytest.approx(result) == 3 / 4 * np.pi

def test_azimuth_third_quarter():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': -100, 'y': -100}
    result = azimuth(first_point, second_point)
    assert pytest.approx(result) == 5 / 4 * np.pi

def test_azimuth_fourth_quarter():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': 100, 'y': -100}
    result = azimuth(first_point, second_point)
    assert pytest.approx(result) == 7 / 4 * np.pi

def test_azimuth_0():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': 100, 'y': 0}
    result = azimuth(first_point, second_point)
    assert pytest.approx(result) == 0

def test_azimuth_100():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': 0, 'y': 100}
    result = azimuth(first_point, second_point)
    assert pytest.approx(result) == np.pi / 2

def test_azimuth_200():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': -100, 'y': 0}
    result = azimuth(first_point, second_point)
    assert pytest.approx(result) == np.pi

def test_azimuth_300():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': 0, 'y': -100}
    result = azimuth(first_point, second_point)
    assert pytest.approx(result) == 3 / 2 * np.pi

def test_azimuth_overlapping_points():
    first_point = {'x': 0, 'y': 0}
    second_point = {'x': 0, 'y': 0}
    with pytest.raises(ValueError):
        azimuth(first_point, second_point)

def test_to_rad_from_degrees():
    result = to_rad(180, 'deg')
    assert pytest.approx(result) == np.pi

def test_to_rad_from_grad():
    result = to_rad(200, 'grad')
    assert pytest.approx(result) == np.pi

def test_to_rad_invalid_unit():
    with pytest.raises(ValueError):
        to_rad(180, 'invalid_unit')

def test_from_rad_to_degrees():
    result = from_rad(np.pi, 'deg')
    assert pytest.approx(result) == 180

def test_from_rad_to_grad():
    result = from_rad(np.pi, 'grad')
    assert pytest.approx(result) == 200

def test_from_rad_invalid_unit():
    with pytest.raises(ValueError):
        from_rad(np.pi, 'invalid_unit')
