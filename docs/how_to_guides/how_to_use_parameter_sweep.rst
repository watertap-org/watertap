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
Depending on how the functions you've defined work, this could be as straightforward as

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

where ``m`` is the flowsheet model that results after the initial "build" step and subsequent operations are performed on that object.

Once this sequence of setup steps is performed, the parameters to be varied should be identified with a dictionary:

.. testcode::

    sweep_params = dict()
    sweep_params['Feed Mass NaCl'] = LinearSample(m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'], 0.005, 0.155, 4)
    sweep_params['Water Recovery'] = LinearSample(m.fs.RO.recovery_mass_phase_comp[0, 'Liq', 'H2O'], 0.3, 0.7, 4)

where the basic pattern is ``dict_name['Short/Pretty-print Name'] = LinearSample(m.path.to.model.variable, lower_limit, upper_limit, num_samples)``.
For example, "Feed Mass NaCl" (the feed mass flow rate of NaCl), which is accessed through the model variable ``m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl']``, is to be varied between 0.005 and 0.155 with 4 equally-spaced values, i.e., ``[0.005, 0.055, 0.105, 0.155]``.
It is also possible to perform random sampling (uniform or normal) with the parameter sweep tool, or the user can specify their own sampling method.
In this case, the 2 parameters will each be varied across 4 values for a total of 16 combinations.
Note that there is no limit on the number of sweep variables specified or their resolution besides the practical limit of how long it will take to optimize using each combination of parameters (e.g., if 5 different variables are provided and each one is individually represented with 20 discrete values, the total number of combinations is 20^5 = 3.2 million!).

After specifying the input parameters, the user should then specify output values on the flowsheet that will be reported in the summary CSV file, which has a similar format to the sweep parameters.
For this RO flowsheet we'll report the levelized cost of water, the optimized RO area, and the output pressure of pump 1:

.. testcode::

    outputs = dict()
    outputs['RO membrane area'] = m.fs.RO.area
    outputs['Pump 1 pressure'] = m.fs.P1.control_volume.properties_out[0].pressure
    outputs['Levelized Cost of Water'] = m.fs.costing.LCOW

Once the problem is setup and the parameters are identified, the parameter_sweep function can finally be invoked which will perform the adjustment and optimization of the model using each combination of variables specified above and saving to `outputs_results.csv` (utilizing the solve method defined in our flowsheet module).
If specified, the full results from each run (the value of every variable and expression) will be reported in `full_results.h5`, along with companion text file containing the metadata of the h5 file in `full_results.txt`.

.. testcode::

    parameter_sweep(m, sweep_params, outputs, csv_results_file='outputs_results.csv', h5_results_file='full_results.h5')

.. testcleanup::

    import os
    os.remove('outputs_results.csv')
    os.remove('full_results.h5')
    os.remove('full_results.txt')

Note that there are additional keyword arguments that can be passed to this function if you desire more control or debugging outputs, especially with regard to the restart logic used after a previous optimization attempt has failed or with managing local outputs computed on parallel hardware.  For more information, consult the technical reference for the parameter sweep tool.

Function Documentation
----------------------

.. automodule :: watertap.tools.parameter_sweep
   :noindex:
   :members:
