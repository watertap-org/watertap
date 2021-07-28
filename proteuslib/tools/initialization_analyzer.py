###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

import weakref
import pyomo.environ as pyo

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from .parameter_sweep import parameter_sweep as _parameter_sweep


class _FlowsheetInitializationAnalyzer:

    def __init__(self, flowsheet, initializer, solver=None):
        self.model = flowsheet
        self._initializer = initializer
        self._store_initial_values()

        if hasattr(self.model, '_fsa'):
            raise RuntimeError(f"Constructed model has attribute _fsa")
        self.model._fsa = weakref.ref(self)
        self.simulate_parameter_sweep = lambda m : m._fsa().simulate()

        if solver is None:
            self._solver = get_solver()
        elif isinstance(solver, str):
            self._solver = pyo.SolverFactory(solver)
        else:
            self._solver = solver

    def simulate(self):
        self._reset_to_initial_values()
        self._initializer(self.model)
        results = self._solver.solve(self.model)
        pyo.assert_optimal_termination(results)

    def _store_initial_values(self):
        m = self.model
        self.initial_value_store = pyo.ComponentMap()
        for var in m.component_data_objects(ctype=pyo.Var, descend_into=True):
            self.initial_value_store[var] = var.value

    def _reset_to_initial_values(self):
        for var, value in self.initial_value_store.items():
            if not var.fixed:
                var.value = value


def initialization_analyzer(model, model_initializer, sweep_params, outputs=None, results_file=None,
        solver=None, mpi_comm=None, num_samples=None, seed=None):
    '''
    This function offers a general way to perform repeated optimizations
    of a model for the purposes of exploring a parameter space while
    monitoring multiple outputs. 
    If provided, writes single CSV file to ``results_file`` with all inputs and resulting outputs.

    Arguments:

        model : A Pyomo ConcreteModel containing a proteuslib flowsheet, for best 
                results it should be initialized before being passed to this
                function.

        model_initializer : A function which takes the provided model as its first argument
                            and performs an initialization routine

        sweep_params: A dictionary containing the values to vary with the format
                      ``sweep_params['Short/Pretty-print Name'] = Sample(model.fs.variable_or_param[index],...)``
                      Samples will be taken as defined in the sampling class.

        outputs (optional) : A dictionary containing "short names" as keys and and Pyomo objects
                             on ``model`` whose values to report as values. E.g.,
                             ``outputs['Short/Pretty-print Name'] = model.fs.variable_or_expression_to_report``.
                             If not provided, a ``feasible`` column will be reported with
                             ``NaN`` if the model was not able to complete initialization or
                             and simulation ``1`` if it was successful

        results_file (optional) : The path and file name where the results are to be saved;
                                   subdirectories will be created as needed.

        mpi_comm (optional) : User-provided MPI communicator for parallel parameter sweeps.
                              If None COMM_WORLD will be used. The default is sufficient for most
                              users.

        num_samples (optional) : If the user is using sampling techniques rather than a linear grid
                                 of values, they need to set the number of samples

        seed (optional) : If the user is using a random sampling technique, this sets the seed

    Returns:

        save_data : A list were the first N columns are the values of the parameters passed
                    by ``sweep_params`` and the remaining columns are the values of the 
                    simulation identified by the ``outputs`` argument.
    '''

    if outputs is None:
        outputs = {'feasible':1}

    fia = _FlowsheetInitializationAnalyzer(model, model_initializer, solver)

    return _parameter_sweep(fia.model, sweep_params, outputs, results_file=results_file, optimize_function=fia.simulate_parameter_sweep,
                            mpi_comm=mpi_comm, num_samples=num_samples, seed=seed)
