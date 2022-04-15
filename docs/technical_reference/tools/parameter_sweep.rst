Parameter Sweep Tool
====================

.. index::
   pair: watertap.tools.parameter_sweep;parameter_sweep

.. currentmodule:: watertap.tools.parameter_sweep

Overview
--------

The parameter sweep tool systematically fixes variables or modifies
mutable parameters on a flowsheet (or Pyomo ConcreteModel), optimizes
the flowsheet, and reports user-specified results.
The parameter sweep tool can be operated in one of two ways: sweeping a
fixed set of parameters, or allowing for random samples from a distribution.
While different fixed or random sampling types can be combined, in a single
parameter sweep the user must use either all fixed or all random sampling types.
(This implementation detail may be relaxed in future releases.)

For the "fixed" sampling type, `LinearSample`, the parameter sweep tool will evaluate the
cross product of all the specified parameters, whereas with the "random" sampling
types `UniformSample` and `NormalSample` the parameter sweep tool will evaluate
a fixed number of samples specified in `num_samples`. With either sampling type,
the values of the Pyomo objects (Vars, Expressions, etc.), that the user specifies
in `outputs`.

For each item the user wants to change, they specify a `sweep_params` dictionary.
The keys are "short" names, and the values are one of the included `Sample` objects.
In all cases the `Sample` objects are instantiated with the Pyomo object to be changed,
with additional arguements depending on the sampling type. For example, for the fixed
`LinearSample` the user would also specify a lower limit, an upper limit, and the number
of elements to be sampled for this parameter between the lower limit and upper limit.
Each item should be sampled at least twice to capture the upper and lower
limit. The random `UniformSample` requires a lower limit and upper limit, and the
`NormalSample` requires a mean and standard deviation.

In addition to the parameters to sweep and the values to track for output,
the user must provide an `optimize_function`, which takes the `model` as an
attribute calls an optimization routine to solve it for the updated parameters.
Should the call to `optimize_function` fail, and a `reinitialize_function` is not
specified, the outputs will be reported as `NaN` for that parameter set.

The user can optionally specify a `reinitialize_function` in case any piece
of the `optimize_function` fails -- after the call to `reinitialize_function`
the model should be in a state ready to optimize again. If the `reinitialize_function`
or the second call to the `optimize_function` fail for any reason, the outputs will
be reported as `NaN` for that parameter set.

The parameter sweep tool maintains the state of the flowsheet / Pyomo model between
calls to `optimize_function` to take advantage of initializations provided by
earlier solutions. If this behavior is undesirable, the user should re-initialize
their flowsheet as part of their `optimize_function`.

Finally, the user can specify a `csv_results_file_name` and/or an `h5_results_file_name`,
which will write the outputs to disk in a CSV and/or H5 format, respectively.
In the CSV results
file, each column specifies a fixed parameter or the associated output, and each row
is a single run with the specified parameters and resulting outputs. The H5 file
contains the parameter sweep inputs and the outputs stored in a dictionary-like format.
Additionally, when an H5 file is written, a companion text file is created with the name
`h5_results_file_name` + `".txt"`. This text file contains the metadata of the H5 results
file.

Parallel Usage
--------------

The parameter sweep tool can optionally utilize `mpi4py` to split the parameter
sweep over multiple processors. If `mpi4py` is installed in the user's conda environment,
a script utilizing the parameter sweep tool can be run in parallel, in this example
using two threads.

.. code:: bash

   mpirun -n 2 python parameter_sweep_script.py

For advanced users, the parameter sweep tool can optionally take a MPI communicator
as an argument.

Function Documentation
----------------------

.. automodule:: watertap.tools.parameter_sweep
   :members:
   :noindex:
