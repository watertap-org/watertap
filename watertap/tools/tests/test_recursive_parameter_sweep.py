###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest
import os
import numpy as np
import pyomo.environ as pyo

from pyomo.environ import value

from watertap.tools.parameter_sweep import (_init_mpi,
                                            _build_combinations,
                                            _divide_combinations,
                                            _update_model_values,
                                            _interp_nan_values,
                                            _process_sweep_params,
                                            _write_output_to_h5,
                                            _read_output_h5,
                                            _create_local_output_skeleton,
                                            _create_global_output,
                                            parameter_sweep,
                                            LinearSample,
                                            UniformSample,
                                            NormalSample,
                                            SamplingType)
from watertap.tools.recursive_parameter_sweep import (_aggregate_filtered_input_arr,
                                                      recursive_parameter_sweep)
from watertap.tools.tests.test_parameter_sweep import _get_rank0_path
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (build,
    set_operating_conditions,
    initialize_system,
    solve,
    optimize)

# import logging
# logging.getLogger('pyomo.core').setLevel(logging.ERROR)

# -----------------------------------------------------------------------------


@pytest.fixture
def model():

    # Initialize a model and solver
    m = pyo.ConcreteModel()

    m.fs = pyo.Block()

    # Declare decision variable and param
    m.fs.x = pyo.Var()
    m.fs.a = pyo.Param(mutable=True)
    m.fs.success_prob = pyo.Param(initialize=0.7)

    # Define expressions and constraints:
    # Numbers must sum to success_prob and x must be positive
    m.fs.err = pyo.Expression(expr = m.fs.a + m.fs.x - m.fs.success_prob)
    m.fs.sum = pyo.Constraint(expr = m.fs.err == 0.0)
    m.fs.pos = pyo.Constraint(expr = m.fs.x >= 0.0)

    return m

    '''
    # Example Usage:
    # Set the number of trials
    nn = 50
    results = np.zeros((nn, 2))

    for k in range(nn):
        # Set a random value for the parameter a
        m.fs.a.set_value(np.random.fs.rand())

        # Attempt to solve the model (infeasible if m.fs.a > success_prob)
        solver.solve(m)

        # Store the value of a and the solution x
        a = pyo.value(m.fs.a)
        x = pyo.value(m.fs.x)

        if np.abs(pyo.value(m.fs.err)) < 1e-6:
            results[k, :] = [a, x]
        else:
            results[k, :] = [a, np.nan]
    '''

@pytest.mark.component
def test_aggregate_filtered_input_arr():
    comm, rank, num_procs = _init_mpi()

    input_dict = {'outputs': {'x_val': {'lower bound': None,
                                        'units': 'None',
                                        'upper bound': None,
                                        'value': np.array([0.31655848, 0.2763452 , 0.26241279, 0.15511682, 0.1511865 , 0.09723662, 0.05410589, 0.6797816 , 0.62896394, 0.6128707 , 0.23852064, 0.17110508, 0.13195544])}},
                  'solve_successful': [True]*13,
                  'sweep_params': {'fs.a': {'units': 'None',
                                            'value': np.array([0.38344152, 0.4236548 , 0.43758721, 0.54488318, 0.5488135 , 0.60276338, 0.64589411, 0.0202184 , 0.07103606, 0.0871293 , 0.46147936, 0.52889492, 0.56804456])}}}


    truth_arr = np.array([0.38344152, 0.4236548 , 0.43758721, 0.54488318, 0.5488135 , 0.60276338, 0.64589411, 0.0202184 , 0.07103606, 0.0871293 , 0.46147936, 0.52889492, 0.56804456])

    req_num_samples = 10
    values_arr = _aggregate_filtered_input_arr(input_dict, req_num_samples, comm, rank, num_procs)
    assert np.allclose(values_arr[:,0], truth_arr[0:req_num_samples], equal_nan=True)

@pytest.mark.component
def test_recursive_parameter_sweep(model, tmp_path):
    comm, rank, num_procs = _init_mpi()
    tmp_path = _get_rank0_path(comm, tmp_path)

    m = model

    solver = pyo.SolverFactory('ipopt')

    sweep_params = {}
    sweep_params['a_val'] = UniformSample(m.fs.a, 0.0, 1.0)

    outputs = {}
    outputs['x_val'] = m.fs.x

    # Run the parameter sweep study using num_samples randomly drawn from the above range
    num_samples = 10
    seed = 0

    results_fname = os.path.join(tmp_path, 'global_results')

    # Run the parameter sweep
    # recursive_parameter_sweep(m, sweep_params, outputs, results_file='recursive_sweep.csv',
    #     optimize_function=optimize, optimize_kwargs={'solver':solver}, req_num_samples=num_samples,
    #     seed=seed, reinitialize_before_sweep=False, reinitialize_function=initialize_system,
    #     reinitialize_kwargs={'solver':solver}
    recursive_parameter_sweep(m, sweep_params, outputs=outputs,
        results_file_name=results_fname, write_csv=False, write_h5=True,
        req_num_samples=num_samples, debugging_data_dir=tmp_path, seed=seed)
