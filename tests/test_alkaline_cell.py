"""This module provides unit tests for `alkaline_cell`"""

import pytest
from numpy.testing import assert_almost_equal

from electrolyzer.simulation.cell_models.alkaline import AlkalineCell as Cell


@pytest.fixture
def cell():
    return Cell.from_dict(
         {
            "model": "default_hri",
            "electrode": {
                "A_electrode": 300,
                "e_e": 0.2,
                "d_em": 0.125,
                "d_ac": 0.25,
            },
            "electrolyte": {
                "w_koh": 30,
            },
            "membrane": {
                "e_m": 0.05,
            },
            "pressure_operating": 1,
            "turndown_ratio": 0.25,
            "max_current_density": 0.3,
            "f_1": 250,
            "f_2": 0.996,
        }
    )

def test_calc_reversible_voltage(cell: Cell):
    """Reversible cell potential should match literature."""
    E_rev = cell.calc_reversible_voltage()

    assert_almost_equal(E_rev, 1.229, decimal=3)