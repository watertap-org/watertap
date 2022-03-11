Monte Carlo testing with the Parameter Sweep
============================================

Overview
--------

This guide introduces two important features of the parameter sweep tool: (1) the ability to prescribe random model values and (2) the ability to take advantage of parallel computing resources for large runs.  It is recommended that new users familiarize themselves with the :ref:`beginner parameter sweep guide<how_to_use_parameter_sweep>` before proceeding.

.. shows you how to use the parameter sweep tool to explore the effect of changing model parameters or decision variables within your WaterTAP model.

.. This might be useful, for example, if you have an existing model of a multi-stage treatment train and you'd like to see the effect of varying Pump 1 pressure and Pump 2 pressure independently (where all possible combinations of Pump 1 and Pump 2 pressure will be explicitly tested).
.. The type and quantity of parameters to be varied are easily changed following steps like the ones below.

How To
------

As before, we begin by importing or explicitly programming any functions relating to flowsheet building/specification, simulation, and optimization setup steps.  We will use the same RO with energy recovery flowsheet for this example.

.. testsetup::

    # quiet idaes logs
    import idaes.logger as idaeslogger
    idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
    idaeslogger.getLogger('ideas.core.util.scaling').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::

    # replace this with your own flowsheet module, e.g.
    # import my_flowsheet_module as mfm
    import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as RO_flowsheet

Once this is done, we import the parameter sweep tool and two different random sampling classes

.. testcode::

    from watertap.tools.parameter_sweep import parameter_sweep, UniformSample, NormalSample

The parameter sweep tool currently offers three random classes:

* ``UniformSample``: Draw samples uniformly from a given upper and lower range.
* ``NormalSample``: Draw samples from a normal distribution given a mean and standard deviation.
* ``LatinHypercubeSample``: Draw samples using a Latin Hypercube algorithm which may yield a more complete exploration of high-dimensional parameter spaces. Note that currently this sample type may not be combined with other sampling types.

We will use the :ref:`same setup steps as before<how_to_use_parameter_sweep>` which returns a flowsheet model ``m``, and performs some initialization

.. testcode::

    # replace these function calls with
    # those in your own flowsheet module

    # set up system
    m = RO_flowsheet.build()
    RO_flowsheet.set_operating_conditions(m)
    RO_flowsheet.initialize_system(m)

    # simulate
    RO_flowsheet.solve(m)

    # set up the model for optimization
    RO_flowsheet.optimize_set_up(m)

.. testoutput::

   ...

Once the model has been setup, we specify the variables to randomly sample using a dictionary

.. testcode::

    sweep_params = dict()
    sweep_params['Spacer_porosity'] = UniformSample(m.fs.RO.spacer_porosity, 0.95, 0.99)
    sweep_params['A_comp'] = NormalSample(m.fs.RO.A_comp, 4.0e-12, 0.5e-12)
    sweep_params['B_comp'] = NormalSample(m.fs.RO.B_comp, 3.5e-8, 0.5e-8)

where the ``spacer_porosity`` attribute will be randomly selected from a uniform distribution of values in the range :math:`[0.95, 0.99]` and model values ``A_comp`` and ``B_comp`` will be drawn from normal distributions centered at :math:`4.0\times10^{-12}` and :math:`3.5\times10^{-8}` with standard deviations of :math:`12-14\%`, respectively.  For this example, we'll extract flowsheet outputs associated with cost, the levelized cost of water (LCOW) and energy consumption (EC), defined via another dictionary

.. testcode::

    outputs = dict()
    outputs['EC'] = m.fs.costing.specific_energy_consumption
    outputs['LCOW'] = m.fs.costing.LCOW


With the flowsheet defined and suitably initialized, along with the definitions for ``sweep_params`` and ``outputs`` on hand, we can call the ``parameter_sweep`` function as before, where we exercise four new keyword arguments: (1) the ability to pass in custom optimization routines to be executed for each sample, (2) the ability to save per-process results for parallel debugging, (3) the specification of the number of samples to draw, and (4) the ability to set a seed for the randomly-generated values which allows consistency to be enforced between runs. The function passed in to `optimize_function` should return a Pyomo results object (i.e., the return value from calling the `solve` method).

.. testcode::

    # Define the local results directory, num_samples, and seed (if desired)
    debugging_data_dir = 'local_results'
    num_samples = 25
    seed = None

    # Run the parameter sweep
    global_results = parameter_sweep(m, sweep_params, outputs, csv_results_file='monte_carlo_results.csv', 
        optimize_function=RO_flowsheet.optimize, debugging_data_dir=debugging_data_dir, num_samples=num_samples, seed=seed)

Note that ``num_samples`` must be provided for any of the random sample classes.  For the very small problem size and simple model used here, parallel hardware is almost certainly not necessary.  However, for larger total numbers of samples or more computationally demanding models, a significant speedup may be attained on a multi-core workstation or high performance computing (HPC) cluster.  To distribute the workload between more than one worker, simply call the scipt using the ``mpirun`` command from the command line

.. code:: bash

    mpirun -n 4 python mc_sweep.py

which will parallelize the requested parameter sweep between 4 computational units, where ``mc_sweep.py`` contains the collection of code snippets shown above ending with the call to ``parameter_sweep``.  Note that there is no requirement that the number of samples be evenly divisible by the number of workers.  In the example shown here with 25 samples and 4 workers, worker 0 processes 7 samples while workers 1-3 process 6 each (you can verify this by examining the four output files in the `local_results` directory).  In most cases, evenly distributing the workload in this way ensures that each worker finishes at roughly the same time.  When each worker has finished, their inidividual results are aggregated into a single result file, `monte_carlo_results.csv`.

.. testcleanup::

    import os
    import shutil
    os.remove('monte_carlo_results.csv')
    shutil.rmtree('local_results')

For more information, consult the technical reference for the parameter sweep tool.

Function Documentation
----------------------

.. automodule :: watertap.tools.parameter_sweep
   :noindex:
   :members:
