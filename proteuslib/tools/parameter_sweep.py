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
import numpy as np
import pyomo.environ as pyo
import sys
import os
import itertools
import warnings

from scipy.interpolate import griddata
from enum import Enum, auto
from abc import abstractmethod, ABC 
from idaes.core.util import get_solver

from idaes.surrogate.pysmo import sampling

# ================================================================

class SamplingType(Enum):
    FIXED = auto()
    RANDOM = auto()
    RANDOM_LHS = auto()

# ================================================================

class _Sample(ABC): 

    def __init__(self, pyomo_object, *args, **kwargs):
        # Check for indexed with single value
        if pyomo_object.is_indexed() and len(pyomo_object) == 1:
            for _data_obj in pyomo_object.values():
                pyomo_object = _data_obj

        # Make sure we are a Var() or Param()
        if not (pyomo_object.is_parameter_type() or pyomo_object.is_variable_type()):
            raise ValueError(f"The sweep parameter needs to be a pyomo Param or Var but {type(pyomo_object)} was provided instead.")
        self.pyomo_object = pyomo_object 
        self.setup(*args, **kwargs)

    @abstractmethod 
    def sample(self, num_samples): 
        pass 

    @abstractmethod 
    def setup(self, *args, **kwargs): 
        pass 

# ================================================================

class RandomSample(_Sample):
    sampling_type = SamplingType.RANDOM

class FixedSample(_Sample):
    sampling_type = SamplingType.FIXED

# ================================================================

class LinearSample(FixedSample):

    def sample(self, num_samples): 
        return np.linspace(self.lower_limit, self.upper_limit, self.num_samples)

    def setup(self, lower_limit, upper_limit, num_samples):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.num_samples = num_samples

# ================================================================

class UniformSample(RandomSample):

    def sample(self, num_samples): 
        return np.random.uniform(self.lower_limit, self.upper_limit, num_samples)

    def setup(self, lower_limit, upper_limit):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit

# ================================================================

class NormalSample(RandomSample):

    def sample(self, num_samples): 
        return np.random.normal(self.mean, self.sd, num_samples)

    def setup(self, mean, sd):
        self.mean = mean
        self.sd = sd

# ================================================================

class LatinHypercubeSample(_Sample):
    sampling_type = SamplingType.RANDOM_LHS

    def sample(self, num_samples): 
        return [self.lower_limit, self.upper_limit]

    def setup(self, lower_limit, upper_limit):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit

# ================================================================

def _init_mpi(mpi_comm=None):

    if mpi_comm is None:
        try:
            from mpi4py import MPI
        except:
            warnings.warn("Could not import mpi4py from current environment (defaulting to serial).")
            return None, 0, 1

        else:
            mpi_comm = MPI.COMM_WORLD

    return mpi_comm, mpi_comm.Get_rank(), mpi_comm.Get_size()

# ================================================================

def _build_combinations(d, sampling_type, num_samples, comm, rank, num_procs):
    num_var_params = len(d)

    if rank == 0:
        param_values = []

        for k, v in d.items():
            # Build a vector of discrete values for this parameter
            p = v.sample(num_samples)
            param_values.append(p)

        if sampling_type == SamplingType.FIXED:
            # Form an array with every possible combination of parameter values
            global_combo_array = np.array(np.meshgrid(*param_values, indexing="ij"))
            global_combo_array = global_combo_array.reshape(num_var_params, -1).T

        elif sampling_type == SamplingType.RANDOM:
            sorting = np.argsort(param_values[0])
            global_combo_array = np.vstack(param_values).T
            global_combo_array = global_combo_array[sorting, :]

        elif sampling_type == SamplingType.RANDOM_LHS:
            lb = [val[0] for val in param_values]
            ub = [val[1] for val in param_values]
            lhs = sampling.LatinHypercubeSampling([lb, ub], number_of_samples=num_samples, sampling_type='creation')
            global_combo_array = lhs.sample_points()
            sorting = np.argsort(global_combo_array[:, 0])
            global_combo_array = global_combo_array[sorting, :]

        else:
            raise ValueError(f"Unknown sampling type: {sampling_type}")

        # Test if the global_combo_array is in row-major order
        if not global_combo_array.flags.c_contiguous:
            # If not, return a copy of this array with row-major memory order
            global_combo_array = np.ascontiguousarray(global_combo_array)

    else:
        if sampling_type == SamplingType.FIXED:
            nx = 1
            for k, v in d.items():
                nx *= v.num_samples
        elif sampling_type == SamplingType.RANDOM or sampling_type == SamplingType.RANDOM_LHS:
            nx = num_samples
        else:
            raise ValueError(f"Unknown sampling type: {sampling_type}")

        if not float(nx).is_integer():
            raise RuntimeError(f"Total number of samples must be integer valued")
        nx = int(nx)

        # Allocate memory to hold the Bcast array
        global_combo_array = np.zeros((nx, num_var_params), dtype=np.float64)

    ### Broadcast the array to all processes
    if num_procs > 1:
        comm.Bcast(global_combo_array, root=0)

    return global_combo_array

# ================================================================

def _divide_combinations(global_combo_array, rank, num_procs):

    # Split the total list of combinations into NUM_PROCS chunks,
    # one per each of the MPI ranks
    # divided_combo_array = np.array_split(global_combo_array, num_procs, axis=0)
    divided_combo_array = np.array_split(global_combo_array, num_procs)

    # Return only this rank's portion of the total workload
    local_combo_array = divided_combo_array[rank]

    return local_combo_array

# ================================================================

def _update_model_values(m, param_dict, values):

    for k, item in enumerate(param_dict.values()):

        param = item.pyomo_object

        if param.is_variable_type():
            # Fix the single value to values[k]
            param.fix(values[k])

        elif param.is_parameter_type():
            # Fix the single value to values[k]
            param.set_value(values[k])

        else:
            raise RuntimeError(f"Unrecognized Pyomo object {param}")

# ================================================================

def _aggregate_results(local_results, global_values, comm, num_procs):

    if num_procs > 1:
        local_results = local_results.astype(np.float64)

        global_results = np.zeros((np.shape(global_values)[0], np.shape(local_results)[1]), dtype=np.float64)

        # Collect the number of result values to be sent from each process
        send_counts = np.zeros(num_procs, dtype=np.int64)
        comm.Gather(np.int64(np.size(local_results)), send_counts, root=0)

        # Collect the global results results onto rank 0
        comm.Gatherv(local_results, (global_results, send_counts), root=0)

        # Broadcast the results to all ranks
        comm.Bcast(global_results, root=0)

    else:
        global_results = np.copy(local_results)

    return global_results

# ================================================================

def _default_optimize(model, options=None, tee=False):
    '''
    Default optimization function used in parameter_sweep.
    Optimizes ``model`` using the IDAES default solver.
    Raises a RuntimeError if the TerminationCondition is not optimal

    Arguments:

        model : A Pyomo ConcreteModel to optimize

        options (optional) : Solver options to pass into idaes.core.utils.get_solver.
                             Default is None
        tee (options) : To display the solver log. Default it False

    '''
    solver = get_solver(options=options)
    results = solver.solve(m, tee=tee)

    if results.solver.termination_condition != pyo.TerminationCondition.optimal:
        raise RuntimeError("The solver failed to converge to an optimal solution. "
                           "This suggests that the user provided infeasible inputs "
                           "or that the model is poorly scaled.")

# ================================================================

def _process_sweep_params(sweep_params):

    sampling_type = None

    # Check the list of parameters to make sure they are valid
    for k in sweep_params:

        # Convert to using Sample class
        if isinstance(sweep_params[k], (list, tuple)):
            sweep_params[k] = LinearSample(*sweep_params[k])

        # Get the type of sampling
        current_sampling_type = sweep_params[k].sampling_type

        # Check to make sure only one sampling type is provided
        if sampling_type is None:
            sampling_type = current_sampling_type
        elif current_sampling_type != sampling_type:
            raise ValueError("Cannot mix sampling types")

    return sweep_params, sampling_type

# ================================================================

def _interp_nan_values(global_values, global_results):

    global_results_clean = np.copy(global_results)

    n_vals = np.shape(global_values)[1]
    n_outs = np.shape(global_results)[1]

    # Build a mask of all the non-nan saved outputs
    # i.e., where the optimzation succeeded
    mask = np.isfinite(global_results[:, 0])

    # Create a list of points where good data is available
    x0 = global_values[mask, :]

    # Interpolate to get a value for nan points where possible
    for k in range(n_outs):
        y0 = global_results[mask, k]
        yi = griddata(x0, y0, global_values, method='linear', rescale=True).reshape(-1)
        global_results_clean[~mask, k] = yi[~mask]

    return global_results_clean

# ================================================================

def parameter_sweep(model, sweep_params, outputs, results_file=None, optimize_function=_default_optimize,
        optimize_kwargs=None, reinitialize_function=None, reinitialize_kwargs=None,
        mpi_comm=None, debugging_data_dir=None, interpolate_nan_outputs=False, num_samples=None, seed=None):

    '''
    This function offers a general way to perform repeated optimizations
    of a model for the purposes of exploring a parameter space while
    monitoring multiple outputs. 
    If provided, writes single CSV file to ``results_file`` with all inputs and resulting outputs.

    Arguments:

        model : A Pyomo ConcreteModel containing a proteuslib flowsheet, for best 
                results it should be initialized before being passed to this
                function.

        sweep_params: A dictionary containing the values to vary with the format
                      ``sweep_params['Short/Pretty-print Name'] =
                      (model.fs.variable_or_param[index], lower_limit, upper_limit, num_samples)``.
                      A uniform number of samples ``num_samples`` will be take between
                      the ``lower_limit`` and ``upper_limit``.

        outputs : A dictionary containing "short names" as keys and and Pyomo objects
                  on ``model`` whose values to report as values. E.g.,
                  ``outputs['Short/Pretty-print Name'] = model.fs.variable_or_expression_to_report``.

        results_file (optional) : The path and file name where the results are to be saved;
                                   subdirectories will be created as needed.

        optimize_function (optional) : A user-defined function to perform the optimization of flowsheet
                                       ``model`` and loads the results back into ``model``. The first
                                       argument of this function is ``model``\. The default uses the
                                       default IDAES solver, raising an exception if the termination
                                       condition is not optimal.

        optimize_kwargs (optional) : Dictionary of kwargs to pass into every call to
                                     ``optimize_function``. The first arg will always be ``model``,
                                     e.g., ``optimize_function(model, **optimize_kwargs)``. The default
                                     uses no kwargs.

        reinitialize_function (optional) : A user-defined function to perform the re-initialize the 
                                           flowsheet ``model`` if the first call to ``optimize_function``
                                           fails for any reason. After ``reinitialize_function``, the
                                           parameter sweep tool will immediately call
                                           ``optimize_function`` again.

        reinitialize_kwargs (optional) : Dictionary or kwargs to pass into every call to 
                                         ``reinitialize_function``. The first arg will always be
                                         ``model``, e.g.,
                                         ``reinitialize_function(model, **reinitialize_kwargs)``.
                                         The default uses no kwargs.

        mpi_comm (optional) : User-provided MPI communicator for parallel parameter sweeps.
                              If None COMM_WORLD will be used. The default is sufficient for most
                              users.

        debugging_data_dir (optional) : Save results on a per-process basis for parallel debugging
                                        purposes. If None no `debugging` data will be saved.

        interpolate_nan_outputs (optional) : When the parameter sweep has finished, interior values
                                             of np.nan will be replaced with a value obtained via
                                             a linear interpolation of their surrounding valid neighbors.
                                             If true, a second output file with the extension "_clean"
                                             will be saved alongside the raw (un-interpolated) values. 

        num_samples (optional) : If the user is using sampling techniques rather than a linear grid
                                 of values, they need to set the number of samples

        seed (optional) : If the user is using a random sampling technique, this sets the seed

    Returns:

        save_data : A list were the first N columns are the values of the parameters passed
                    by ``sweep_params`` and the remaining columns are the values of the 
                    simulation identified by the ``outputs`` argument.
    '''

    # Get an MPI communicator
    comm, rank, num_procs = _init_mpi(mpi_comm)

    # Convert sweep_params to LinearSamples
    sweep_params, sampling_type = _process_sweep_params(sweep_params)

    # Set the seed before sampling 
    np.random.seed(seed)

    # Enumerate/Sample the parameter space
    global_values = _build_combinations(sweep_params, sampling_type, num_samples, comm, rank, num_procs)

    # divide the workload between processors
    local_values = _divide_combinations(global_values, rank, num_procs)

    # Initialize space to hold results
    local_num_cases = np.shape(local_values)[0]
    local_results = np.zeros((local_num_cases, len(outputs)))

    # Set up optimize_kwargs
    if optimize_kwargs is None:
        optimize_kwargs = dict()
    # Set up reinitialize_kwargs
    if reinitialize_kwargs is None:
        reinitialize_kwargs = dict()

    # ================================================================
    # Run all optimization cases
    # ================================================================

    for k in range(local_num_cases):
        # Update the model values with a single combination from the parameter space
        _update_model_values(model, sweep_params, local_values[k, :])

        try:
            # Simulate/optimize with this set of parameters
            optimize_function(model, **optimize_kwargs)

        except:
            # If the run is infeasible, report nan
            local_results[k, :] = np.nan
            previous_run_failed = True

        else:
            # If the simulation suceeds, report stats
            local_results[k, :] = [pyo.value(outcome) for outcome in outputs.values()]
            previous_run_failed = False

        if previous_run_failed and (reinitialize_function is not None):
            # We choose to re-initialize the model at this point
            try:
                reinitialize_function(model, **reinitialize_kwargs)
                optimize_function(model, **optimize_kwargs)
            except:
                # do we raise an error here?
                # nothing to do
                pass
            else:
                local_results[k, :] = [pyo.value(outcome) for outcome in outputs.values()]


    # ================================================================
    # Save results
    # ================================================================

    global_results = _aggregate_results(local_results, global_values, comm, num_procs)

    # Make a directory for saved outputs
    if rank == 0:
        if results_file is not None:
            dirname = os.path.dirname(results_file)
            if dirname != '':
                os.makedirs(dirname, exist_ok=True)
                
        if debugging_data_dir is not None:
            os.makedirs(debugging_data_dir, exist_ok=True)

    if num_procs > 1:
        comm.Barrier()

    # Write a header string for all data files
    data_header = ','.join(itertools.chain(sweep_params,outputs))

    if debugging_data_dir is not None:
        # Create the local filename and data
        fname = os.path.join(debugging_data_dir, f'local_results_{rank:03}.csv')
        local_save_data = np.hstack((local_values, local_results))

        # Save the local data
        np.savetxt(fname, local_save_data, header=data_header, delimiter=', ', fmt='%.6e')

    # Create the global filename and data
    global_save_data = np.hstack((global_values, global_results))

    if rank == 0 and results_file is not None:
        # Save the global data
        np.savetxt(results_file, global_save_data, header=data_header, delimiter=',', fmt='%.6e')

        if interpolate_nan_outputs:
            global_results_clean = _interp_nan_values(global_values, global_results)
            global_save_data_clean = np.hstack((global_values, global_results_clean))

            head, tail = os.path.split(results_file)

            if head == '':
                interp_file = 'interpolated_%s' % (tail)
            else:
                interp_file = '%s/interpolated_%s' % (head, tail)

            np.savetxt(interp_file, global_save_data_clean, header=data_header, delimiter=',', fmt='%.6e')
    
    return global_save_data

# ================================================================
