How to Explore a Model with Parameter Sweep
===========================================

Overview
--------

This guide shows you how to use the parameter sweep tool to explore the effect of changing model parameters or decision variables within your ProteusLib model.  This is might be useful, for example, if you have an existing model of a multi-stage treatment train and you'd like to see the effect of varying Pump 1 pressure and Pump 2 pressure independently (where all possible combinations of Pump 1 and Pump 2 pressure will be explicitly tested).  The type and quantity of parameters to be varied are easily changed following steps like the ones below.

How To
------

Begin by importing or explicitly programming any functions relating to the following steps

1. Building the Pyomo flowsheet model
2. Simulating the model using an initial condition
3. Setting up the optimization (e.g., unfixing and setting bounds)
4. Performing the optimization

For example, the code below imports the 4 critical functions from an existing flowsheet ::

    from LSRRO_model import build_model, simulate, set_up_optimization, optimize

Once this is done, import the parameter sweep tool ::

    from proteuslib.tools.parameter_sweep import parameter_sweep

Conceptually, regardless of the number of iterations necessary to test each possible combination of variables, it is only necessary to build, simulate, and set up the model once.  Thus, these steps are left to the user and handled outside the parameter sweep function.  Depending on how the functions you've defined work, this could be as straightforward as ::

    m = build_model(optional_args=...)
    m = simulate(m, optional_args=...)
    m = set_up_optimization(m, optional_args=...)

where ``m`` is the flowsheet model that results after the initial "build" step and subsequent operations are performed on that object.  Once this sequence of setup steps is performed, the parameters to be varied should be identified with a dictionary. ::

    sweep_params = dict()
    sweep_params['Feed Concentration'] = (m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'], 0.005, 0.155, 4)
    sweep_params['Water Recovery'] = (m.fs.water_recovery, 0.3, 0.7, 4)
    ...

where the basic pattern is ``dict_name['Short/Pretty-print Name'] = (m.path.to.model.variable, lower_limit, upper_limit, num_samples)``.  For example, "Feed Concentration", which is accessed in the model variable ``m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']``, is to be varied between 0.005 and 0.155 with 4 equally-spaced values, i.e., ``[0.005, 0.055, 0.105, 0.155]``.  In this case, the 2 parameters will each be varied across 4 values for a total of 16 combinations.  Note that there is no limit on the number of sweep variables specified or their resolution besides the practical limit of how long it will take to optimize using each combination of parameters (e.g., if 5 different variables are provided and each one is individually represented with 20 discrete values, the total number of combinations is 20^5 = 3.2 million!). Once the problem is setup and the parameters are identified, the parameter_sweep function can finally be invoked which will perform the adjustment and optimization of the model invoking each combination of variables specified above. ::

    parameter_sweep(m, sweep_params, outputs, path/to/results/filename.csv, optimize, optimize_kwargs={...}):

Note that there are additional keyword arguments that can be passed to this function if you desire more control or debugging outputs, especially with regard to the restart logic used after a previous optimization attempt has failed or with managing local outputs computed on parallel hardware.  For more information, consult the technical reference for the parameter sweep tool.

Function Documentation
----------------------

.. automodule :: proteuslib.tools.parameter_sweep
   :members:
