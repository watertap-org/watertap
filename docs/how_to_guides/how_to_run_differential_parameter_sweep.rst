How to Run Differential Parameter Sweep
============================================

Overview
--------

This guide explains how to run a differential parameter sweep for conducting, e.g., value of innovation (VOI), analysis.

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

Once this is done, we import the differential parameter sweep tool and sampling classes.

.. testcode::

    from watertap.tools.parameter_sweep import differential_parameter_sweep, UniformSample, NormalSample

We will use the :ref:`same setup steps as before<how_to_use_parameter_sweep>` to set up the generating functions for our model, sweep params, and outputs:

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

Once the model has been setup, we specify the variables to sample using a dictionary

.. testcode::

    def build_sweep_params(model, num_samples=5):
        sweep_params = dict()
        sweep_params["A_comp"] = UniformSample(model.fs.RO.A_comp, 4.2e-12, 2.1e-11, num_samples)
        sweep_params["membrane_cost"] = UniformSample(
            model.fs.costing.reverse_osmosis.membrane_cost, 30, 10, num_samples
        )
        return sweep_params

where the ``A_comp`` attribute will be randomly selected from a uniform distribution of values in the range :math:`[4.2e-12, 2.1e-11]` and ``membrane_cost`` will be drawn from a uniform distribution between :math:`[30, 10]`.  For this example, we will extract flowsheet outputs associated with cost, the levelized cost of water (LCOW) and energy consumption (EC), defined via another dictionary.
Next we will construct a specification dictionary to run the differential parameter sweep.

.. testcode::

    def build_diff_sweep_param_specs(model):
        differential_sweep_specs = dict()
        
        differential_sweep_specs["A_comp"] = {
            "diff_sample_type": NormalSample,
            "std_dev": 0.3e-12,
            "pyomo_object": model.fs.RO.A_comp,
        }

        differential_sweep_specs["membrane_cost"] = {
            "diff_sample_type": UniformSample,
            "diff_mode": "percentile",
            "nominal_lb" : sweep_params["membrane_cost"].lower_limit,
            "nominal_ub" : sweep_params["membrane_cost"].upper_limit,
            "relative_lb" : -0.05,
            "relative_ub" : -0.05,
            "pyomo_object": model.fs.costing.reverse_osmosis.membrane_cost,
        }

        return differential_sweep_specs

``differential_sweep_specs`` is a specification dictionary that contains details for how to construct the parameter sweep dictionary for differential sweep. This is a nested dictionary where the first level denotes the variable names for which the differential sweep needs to be carried out. The second level denotes various options to be used for each variable. The number of samples for each differential sweep is specified while initializing the ``DifferentialParameterSweep`` object using the keyword ``num_diff_samples``. There are 4 modes of setting up a variable to undergo differential sweep:

#. ``NormalSample`` : Uses the nominal value as the mean and expects ``std_dev`` key for the differential sweep sampling. It looks like the following:

    .. code-block:: python

        differential_sweep_specs["A_comp"] = {
                "diff_sample_type": NormalSample,
                "std_dev": 0.3e-12,
                "pyomo_object": model.fs.RO.A_comp,
            }

    This differential mode is unique to variables that expect normal sampling. *All other sampling types expect one of the other 3 differential modes below.*

#. ``sum`` : Perturbs the nominal value by a certain absolute percentage to create an upper and lower bound for the differential solve. The logic in the code looks as follows:

    .. code-block:: python

        lower_bound = nominal_val * (1 - relative_lb)
        upper_bound = nominal_val * (1 + relative_ub)

#. ``product``: Perturbs the nominal value by a scaling factor to create upper and lower bounds for the differential sweep. It uses the following logic

    .. code-block:: python

        lower_bound = nominal_val * relative_lb
        upper_bound = nominal_val * relative_ub  

#. ``percentile``: Perturbs the nominal value by a percentage of the difference between the nominal upper and lower bound values. The logic is 

    .. code-block:: python

        delta_nominal = abs(upper_nominal - lower_nominal)
        lower_bound = nominal_val + delta_nominal * relative_lb
        upper_bound = nominal_val + delta_nominal * relative_ub

An example differential sweep spec dictionary may look like the following:

.. code-block:: python

    differential_sweep_specs = dict()
    differential_sweep_specs["membrane_cost"] = {
            "diff_sample_type": UniformSample,
            "diff_mode": "percentile",
            "nominal_lb" : sweep_params["membrane_cost"].lower_limit,
            "nominal_ub" : sweep_params["membrane_cost"].upper_limit,
            "relative_lb" : -0.05,
            "relative_ub" : -0.05,
            "pyomo_object": model.fs.costing.reverse_osmosis.membrane_cost,
        }
    differential_sweep_specs["px_cost"] = {
        "diff_sample_type": LinearSample,
        "diff_mode": "sum",
        "relative_lb" : -0.05,
        "relative_ub" : -0.05,
        "pyomo_object": m.fs.costing.pressure_exchanger.cost,
    }
    differential_sweep_specs["px_efficiency"] = {
        "diff_sample_type": UniformSample,
        "diff_mode": "product",
        "relative_lb" : 0.001,
        "relative_ub" : 0.001,
        "pyomo_object": m.fs.PXR.efficiency_pressure_exchanger,
    }

.. important:: The user can only conduct differential sweeps for variables specified with ``sweep_params``.

Continuing with the example test code from above, we will use the following function for building the outputs.

.. testcode::

    def build_outputs(model, sweep_params):
        outputs = dict()
        outputs['EC'] = model.fs.costing.specific_energy_consumption
        outputs['LCOW'] = model.fs.costing.LCOW
        return outputs

With the flowsheet defined and suitably initialized, along with the definitions for ``sweep_params``, ``differential_sweep_specs``, and ``outputs`` on hand, we can call the ``differential_parameter_sweep`` function as before.

.. note:: This documentation currently uses the older API for calling the differential parameter sweep. This API will be deprecated in the near future. The documentation will be changed to reflect this accordingly. We recommend running the differential parameter sweep in serial or with MPI only.

.. testcode::

    # Define the local results directory, num_samples, and seed (if desired)
    num_samples = 5
    seed = None

    model = build_model()
    sweep_params = build_sweep_params(model, num_samples=num_samples)
    differential_sweep_specs = build_diff_sweep_param_specs(model)
    outputs = build_outputs(model, sweep_params)

    # Run the parameter sweep
    global_results = differential_parameter_sweep(
            build_model, 
            build_sweep_params, 
            differential_sweep_specs,
            outputs, 
            h5_results_file_name='monte_carlo_results.h5',
            optimize_function=RO_flowsheet.optimize,
            debugging_data_dir=None,
            num_samples=num_samples,
            num_diff_samples=2,
            seed=seed,
        )

.. testoutput::

    ...

.. testcleanup::

    import os
    import shutil
    try:
        os.remove('monte_carlo_results.h5')
        os.remove('monte_carlo_results.h5.txt')
    except:
        print("monte_carlo_results.h5 does not exist, nothing to delete.")

Module Documentation
--------------------

* :mod:`watertap.tools.parameter_sweep`