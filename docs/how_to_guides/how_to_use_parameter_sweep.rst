.. _how_to_use_parameter_sweep:

How to explore a model with parameter sweep
===========================================

Overview
--------

This guide shows you how to use the parameter sweep tool to explore the effect of changing model parameters or decision variables within your WaterTAP model.
This might be useful, for example, if you have an existing model of a multi-stage treatment train and you'd like to see the effect of varying Pump 1 pressure and Pump 2 pressure independently (where all possible combinations of Pump 1 and Pump 2 pressure will be explicitly tested).
The type and quantity of parameters to be varied are easily changed following steps like the ones below.

How To
------

Begin by importing or explicitly programming any functions relating to the following steps:

1. Building the Pyomo flowsheet model
2. Simulating the model using an initial condition
3. Setting up the optimization (e.g., unfixing and setting bounds)
4. Performing the optimization

For example, the code below imports from an existing flowsheet module, RO with energy recovery.
In general you would import your own flowsheet module.

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

Once this is done, import the parameter sweep tool

.. testcode::

    from watertap.tools.parameter_sweep import parameter_sweep, LinearSample

Conceptually, regardless of the number of iterations necessary to test each possible combination of variables, it is only necessary to build, simulate, and set up the model once.
Thus, these steps are left to the user and handled outside the parameter sweep function.

In order to support native parallelism within the parameter sweep class, the preferred way of setting up a model is to create a standalone function 
that can produce it. Depending on how the functions you've defined work, this could be as straightforward as

.. testcode::

    def build_model(**kwargs):

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

        return m

where ``m`` is the flowsheet model that results after the initial "build" step and subsequent operations are performed on that object.

Once this sequence of setup steps is performed, the parameters to be varied should be identified with a dictionary. Similarly to the way a 
model is produced, in order to support native parallelism the preferred way to define parameters is by defining a function that creates them
(the generated model will automatically be passed in as the first argument):

.. testcode::

    def build_sweep_params(model, **kwargs):
        sweep_params = dict()
        sweep_params['Feed Mass NaCl'] = LinearSample(model.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'], 0.005, 0.155, 4)
        sweep_params['Water Recovery'] = LinearSample(model.fs.RO.recovery_mass_phase_comp[0, 'Liq', 'H2O'], 0.3, 0.7, 4)
        return sweep_params

where the basic pattern is ``dict_name['Short/Pretty-print Name'] = LinearSample(m.path.to.model.variable, lower_limit, upper_limit, num_samples)``.
For example, "Feed Mass NaCl" (the feed mass flow rate of NaCl), which is accessed through the model variable ``m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl']``, is to be varied between 0.005 and 0.155 with 4 equally-spaced values, i.e., ``[0.005, 0.055, 0.105, 0.155]``.
It is also possible to perform random sampling (uniform or normal) with the parameter sweep tool, or the user can specify their own sampling method.
In this case, the 2 parameters will each be varied across 4 values for a total of 16 combinations.
Note that there is no limit on the number of sweep variables specified or their resolution besides the practical limit of how long it will take to optimize using each combination of parameters (e.g., if 5 different variables are provided and each one is individually represented with 20 discrete values, the total number of combinations is 20^5 = 3.2 million!).

After specifying the input parameters, the user should then specify output values on the flowsheet that will be reported in the summary CSV file, which has a similar format to the sweep parameters.
For this RO flowsheet we'll report the levelized cost of water, the optimized RO area, and the output pressure of pump 1:

.. testcode::

    def build_outputs(model, **kwargs):
        outputs = dict()
        outputs['RO membrane area'] = model.fs.RO.area
        outputs['Pump 1 pressure'] = model.fs.P1.control_volume.properties_out[0].pressure
        outputs['Levelized Cost of Water'] = model.fs.costing.LCOW
        return outputs

Once the problem is setup and the parameters are identified, the parameter_sweep function can finally be invoked which will perform the adjustment and optimization of the model using each combination of variables specified above (utilizing the solve method defined in our flowsheet module).
If specified, the parameter_sweep function will optionally write results in CSV format to the path specified in `csv_results_file_name` or in H5 format to the path specified in `h5_results_file_name`.
The file `outputs_results.csv` contains the `sweep_param` values and `outputs` values in an array format, while `outputs_results.h5` contains a dictionary containing the `sweep_params`, `outputs`, and a boolean list of successful or failed solves.
The H5 writer also creates a companion text file containing the metadata of the h5 file in `outputs_results.h5.txt`.
Note that if `outputs = None` and an H5 results file is specified, all of the pyomo model variables will be stored in the `outputs_results.h5` and `outputs_results.h5.txt` files.

Passing in a model, sweep params, and outputs directly to the parameter_sweep function is currently supported but is deprecated and will be removed in
future versions. The preferred way is to pass in generating functions as shown below:

.. testcode::

    parameter_sweep(build_model, build_sweep_params, build_outputs, csv_results_file_name='outputs_results.csv', h5_results_file_name='outputs_results.h5')

.. testoutput::

    ...

.. testcleanup::

    import os
    os.remove('outputs_results.csv')
    os.remove('outputs_results.h5')
    os.remove('outputs_results.h5.txt')

Note that there are additional keyword arguments that can be passed to this function if you desire more control or debugging outputs, especially with regard to the restart logic used after a previous optimization attempt has failed or with managing local outputs computed on parallel hardware.  For more information, consult the technical reference for the parameter sweep tool.

Module Documentation
--------------------

* :mod:`watertap.tools.parameter_sweep`
