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

from watertap.tools.parameter_sweep.parameter_sweep import (
    ParameterSweep,
    RecursiveParameterSweep,
)


def parameter_sweep(
    model,
    sweep_params,
    outputs=None,
    csv_results_file_name=None,
    h5_results_file_name=None,
    optimize_function=None,
    optimize_kwargs=None,
    reinitialize_function=None,
    reinitialize_kwargs=None,
    reinitialize_before_sweep=False,
    probe_function=None,
    mpi_comm=None,
    debugging_data_dir=None,
    interpolate_nan_outputs=False,
    num_samples=None,
    seed=None,
):

    """
    This function offers a general way to perform repeated optimizations
    of a model for the purposes of exploring a parameter space while
    monitoring multiple outputs.
    If provided, writes single CSV file to ``results_file`` with all inputs and resulting outputs.

    Arguments:

        model : A Pyomo ConcreteModel containing a watertap flowsheet, for best
                results it should be initialized before being passed to this
                function.

        sweep_params: A dictionary containing the values to vary with the format
                      ``sweep_params['Short/Pretty-print Name'] =
                      (model.fs.variable_or_param[index], lower_limit, upper_limit, num_samples)``.
                      A uniform number of samples ``num_samples`` will be take between
                      the ``lower_limit`` and ``upper_limit``.

        outputs : An optional dictionary containing "short names" as keys and and Pyomo objects
                  on ``model`` whose values to report as values. E.g.,
                  ``outputs['Short/Pretty-print Name'] = model.fs.variable_or_expression_to_report``.
                  If not provided, i.e., outputs = None, the default behavior is to save all model
                  variables, parameters, and expressions which provides very thorough results
                  at the cost of large file sizes.

        csv_results_file_name (optional) : The path and file name to write a csv file. The default `None`
                                           does not write a csv file.

        h5_results_file_name (optional) : The path and file name to write a h5 file. The default `None`
                                          does not write a file.
                                          Writing an h5 file will also create a companion text file `{h5_results_file_name}.txt`
                                          which contains the variable names contained within the H5 file.

        optimize_function (optional) : A user-defined function to perform the optimization of flowsheet
                                       ``model`` and loads the results back into ``model``. The first
                                       argument of this function is ``model``. The default uses the
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

        reinitialize_before_sweep (optional): Boolean option to reinitialize the flow sheet model before
                                              every parameter sweep realization. The default is False.
                                              Note the parameter sweep model will try to reinitialize the
                                              solve regardless of the option if the run fails.

        probe_function (optional): A user-defined function that can cheaply check if a current model
                                  configuration is solvable without actually reinitializing or solving.

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
    """

    kwargs = {}
    if mpi_comm is not None:
        kwargs["comm"] = mpi_comm
    if csv_results_file_name is not None:
        kwargs["csv_results_file_name"] = csv_results_file_name
    if h5_results_file_name is not None:
        kwargs["h5_results_file_name"] = h5_results_file_name
    if optimize_function is not None:
        kwargs["optimize_function"] = optimize_function
    if optimize_kwargs is not None:
        kwargs["optimize_kwargs"] = optimize_kwargs
    if reinitialize_function is not None:
        kwargs["reinitialize_function"] = reinitialize_function
    if reinitialize_kwargs is not None:
        kwargs["reinitialize_kwargs"] = reinitialize_kwargs
    kwargs["reinitialize_before_sweep"] = reinitialize_before_sweep
    if probe_function is not None:
        kwargs["probe_function"] = probe_function
    if debugging_data_dir is not None:
        kwargs["debugging_data_dir"] = debugging_data_dir
    if interpolate_nan_outputs is not None:
        kwargs["interpolate_nan_outputs"] = interpolate_nan_outputs

    ps = ParameterSweep(**kwargs)

    return ps.parameter_sweep(
        model, sweep_params, outputs=outputs, num_samples=num_samples, seed=seed
    )


def recursive_parameter_sweep(
    model,
    sweep_params,
    outputs=None,
    csv_results_file_name=None,
    h5_results_file_name=None,
    optimize_function=None,
    optimize_kwargs=None,
    reinitialize_function=None,
    reinitialize_kwargs=None,
    reinitialize_before_sweep=False,
    probe_function=None,
    mpi_comm=None,
    debugging_data_dir=None,
    interpolate_nan_outputs=False,
    req_num_samples=None,
    seed=None,
):

    kwargs = {}
    if mpi_comm is not None:
        kwargs["comm"] = mpi_comm
    if csv_results_file_name is not None:
        kwargs["csv_results_file_name"] = csv_results_file_name
    if h5_results_file_name is not None:
        kwargs["h5_results_file_name"] = h5_results_file_name
    if optimize_function is not None:
        kwargs["optimize_function"] = optimize_function
    if optimize_kwargs is not None:
        kwargs["optimize_kwargs"] = optimize_kwargs
    if reinitialize_function is not None:
        kwargs["reinitialize_function"] = reinitialize_function
    if reinitialize_kwargs is not None:
        kwargs["reinitialize_kwargs"] = reinitialize_kwargs
    kwargs["reinitialize_before_sweep"] = reinitialize_before_sweep
    if probe_function is not None:
        kwargs["probe_function"] = probe_function
    if debugging_data_dir is not None:
        kwargs["debugging_data_dir"] = debugging_data_dir
    if interpolate_nan_outputs is not None:
        kwargs["interpolate_nan_outputs"] = interpolate_nan_outputs

    rps = RecursiveParameterSweep(**kwargs)

    return rps.parameter_sweep(
        model, sweep_params, outputs=outputs, req_num_samples=req_num_samples, seed=seed
    )
