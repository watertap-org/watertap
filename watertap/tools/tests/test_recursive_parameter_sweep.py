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
                                            # _aggregate_results,
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
from watertap.tools.recursive_parameter_sweep import *
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

    # Run the parameter sweep
    # recursive_parameter_sweep(m, sweep_params, outputs, results_file='recursive_sweep.csv',
    #     optimize_function=optimize, optimize_kwargs={'solver':solver}, req_num_samples=num_samples,
    #     seed=seed, reinitialize_before_sweep=False, reinitialize_function=initialize_system,
    #     reinitialize_kwargs={'solver':solver}
    recursive_parameter_sweep(m, sweep_params, outputs=outputs,
        results_dir=tmp_path, results_fname='recursive_output',
        req_num_samples=num_samples, seed=seed)
