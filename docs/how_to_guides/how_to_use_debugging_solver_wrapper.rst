.. _how_to_use_debugging solver wrapper:

How to use the debugging solver wrapper
=======================================

Overview
--------

Debugging a failed initialization or solve can be cumbersome. 
Further, it is often useful to have the model state (i.e., initial variable values) from before the solve failed.
The debugging solver wrapper facilitates the following:

(1) Stores initialization information before a failed solve

(2) Upon a failed attempt to solve, the user is routed to an Interactive Python notebook where the restored model state can be accessed, `IDAES' DiagnosticToolbox <https://idaes-pse.readthedocs.io/en/stable/reference_guides/core/util/diagnostics/diagnostics_toolbox.html>`_ is instantiated to probe the model, and the user can freely apply any other diagnostic utility functions to troubleshoot the problematic model. 

How To
------

In a python module containing the model and script to solve that model, the user would make a simple import:

.. testcode::

    from watertap_solvers.model_debug_mode import activate
    activate()


.. warning::
 
    If you ran your python file in an interactive window, this debugging mode may not work as expected. We recommend running your python file in a terminal.

Example behavior without debugging mode
---------------------------------------

The example output below shows a problematic model that fails to initialize.

.. code-block:: text

    2024-02-01 13:06:07 [DEBUG] idaes.solve.fs.bed_stack: EXIT: Converged to a point of local infeasibility. Problem may be infeasible.
    2024-02-01 13:06:07 [DEBUG] idaes.solve.fs.bed_stack: WARNING: Loading a SolverResults object with a warning status into
    2024-02-01 13:06:07 [DEBUG] idaes.solve.fs.bed_stack: model.name="fs.bed_stack";
    2024-02-01 13:06:07 [DEBUG] idaes.solve.fs.bed_stack:     - termination condition: infeasible
    2024-02-01 13:06:07 [DEBUG] idaes.solve.fs.bed_stack:     - message from solver: Ipopt 3.14.11\x3a Converged to a locally infeasible
    2024-02-01 13:06:07 [DEBUG] idaes.solve.fs.bed_stack:       point. Problem may be infeasible.
    2024-02-01 13:06:07 [INFO] idaes.init.fs.bed_stack: Initialization Step 3 infeasible - Converged to a locally infeasible point. Problem may be infeasible..
    2024-02-01 13:06:07 [WARNING] idaes.init.fs.bed_stack:  The solver at the Initialization Step 3 step failed to converge to an optimal solution.This suggests that the user provided infeasible inputs or that the model is poorly scaled, poorly initialized, or degenerate.
    2024-02-01 13:06:07 [INFO] idaes.init.fs.bed_stack: Initialization Complete: infeasible - Converged to a locally infeasible point. Problem may be infeasible.
    Traceback (most recent call last):
    File "/Models/bed_simulation.py", line 439, in <module>
    m, res = main()
    File "/Models/bed_simulation.py", line 53, in main
    m, res = run_simulation(case, parameter_estimates)
    File "/Models/bed_simulation.py", line 106, in run_simulation
    model_initialize(m, case)
    File "/Models/bed_simulation.py", line 313, in model_initialize
    model.fs.bed_stack.initialize(outlvl=idaeslog.DEBUG, ignore_dof=True)
    File "/watertap/core/initialization_mixin.py", line 23, in initialize
    return super().initialize(*args, **kwargs)
    File "/anaconda3/envs/watertap/lib/python3.10/site-packages/idaes/core/base/unit_model.py", line 540, in initialize
    flags = blk.initialize_build(*args, **kwargs)
    File "/watertap/unit_models/electrodialysis_1D.py", line 2146, in initialize_build
    raise InitializationError(f"Unit model {self.name} failed to initialize")
    idaes.core.util.exceptions.InitializationError: Unit model fs.bed_stack failed to initialize

Example behavior with debugging mode
---------------------------------------
Adding the aforementioned import to the module and calling ``activate()`` results in the printout below before being routed to an Interactive Python window:

.. code-block:: text

    EXIT: Converged to a point of local infeasibility. Problem may be infeasible.
    WARNING: Loading a SolverResults object with a warning status into
    model.name="fs.bed_stack";
     - termination condition: infeasible
     - message from solver: Ipopt 3.14.11\x3a Converged to a locally infeasible
        point. Problem may be infeasible.
    
    Solver debugging mode: the block fs.bed_stack failed to solve.
    
    fs.bed_stack can be called as `blk` in debugging mode.
    
    The solver ipopt-watertap is available in the variable `solver`.
    
    The initial values before the failed solve have been stored.
    
    You can restore these initial values at anytime by calling `debug.restore_initial_values(blk)`.
    
    The model has been loaded into an IDAES DiagnosticsToolbox instance called `dt`.
    
    WARNING: If you ran your python file in an interactive window, this debugging mode will not work as intended. Be sure to run your python file in a terminal.
    
    Python 3.10.9 (main, Jan 11 2023, 09:18:20) [Clang 14.0.6 ]
    Type 'copyright', 'credits' or 'license' for more information
    IPython 7.34.0 -- An enhanced Interactive Python. Type '?' for help.

    In [1]: 

Check the model name with ``blk``:

.. code-block:: shell

    In [1]: blk.name
    Out[1]: 'fs.bed_stack'

Use the DiagnosticsToolbox (instantiated to ``dt``) to probe for structural issues in the model:

.. code-block:: shell

    In [2]: dt.report_structural_issues()
    ====================================================================================
    Model Statistics

            Activated Blocks: 15 (Deactivated: 0)
            Free Variables in Activated Constraints: 566 (External: 0)
                Free Variables with only lower bounds: 136
                Free Variables with only upper bounds: 0
                Free Variables with upper and lower bounds: 240
            Fixed Variables in Activated Constraints: 42 (External: 7)
            Activated Equality Constraints: 566 (Deactivated: 0)
            Activated Inequality Constraints: 0 (Deactivated: 0)
            Activated Objectives: 0 (Deactivated: 0)

    ------------------------------------------------------------------------------------
    1 WARNINGS

        WARNING: Found 354 potential evaluation errors.

    ------------------------------------------------------------------------------------
    2 Cautions

        Caution: 3 variables fixed to 0
        Caution: 11 unused variables (0 fixed)

    ------------------------------------------------------------------------------------
    Suggested next steps:

        display_potential_evaluation_errors()

    ====================================================================================

Continue to probe and diagnose model infeasibility in this Interactive Python window.
