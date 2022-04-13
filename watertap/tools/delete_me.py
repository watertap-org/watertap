Current:

sweep_params = {} # Dictionary with 3 parameters
n_sweep_params = 3
n_samples = 100
initial_input = np.zeros((3,100))

while no_samples:
    for i in range(num_samples):
        result = solve(initial_input[i,:])

Proposed:

n_sweep_params = 3
differential_mc_idx = 1 # We will know this
n_samples = 100
initial_input = np.zeros((3,100))

# Changes will be by percentile or functions

for i in range(num_samples):
    result = solve(initial_input[i,:])
    # if result == some_condition:
    diff_samples = sample_diff_func()
    n_sub_samples = len(diff_samples)
    for j in range(n_sub_samples):
        sub_result = solve(diff_samples[j])

# Two different loops that loads in the results.
# Hot restart?



################


_ParameterSweepBase
ParameterSweep <: _ParameterSweepBase
RecursivePrameterSweep <: _ParameterSweepBase
DifferentialPrameterSweep <: _ParameterSweepBase

_ParameterSweepWriter
_ParameterSweepReader

"""
Using `pyomo_config_objects` for options. Anyone iusing idaes will be very
familiar
https://pyomo.readthedocs.io/en/stable/developer_reference/config.html
"""
