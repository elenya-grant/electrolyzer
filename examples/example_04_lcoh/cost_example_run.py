"""
This example performs a fully controlled electrolyzer simulation using the
`run_electrolyzer` function. See `example_run.ipynb` for an interactive option.
"""
import os

import numpy as np

from electrolyzer import run_lcoh  # , run_electrolyzer,LCOH,


fname_input_modeling = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "cost_modeling_options.yaml"
)

turbine_rating = 3.4  # MW

# Create cosine test signal
test_signal_angle = np.linspace(0, 8 * np.pi, 3600 * 8 + 10)
base_value = (turbine_rating / 2) + 0.2
variation_value = turbine_rating - base_value
power_test_signal = (base_value + variation_value * np.cos(test_signal_angle)) * 1e6

# res = run_electrolyzer(fname_input_modeling, power_test_signal)
lcoe = 44.18 * (1 / 1000)
res = run_lcoh(fname_input_modeling, power_test_signal, lcoe)
#
# run_lcoh(fname_input_modeling, res[0], res[1], lcoe)
[]
# print(res)
