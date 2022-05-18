WaterTAP User Interface
========================

.. py:currentmodule:: watertap.ui.api

Overview
--------

The WaterTAP user interface is an application programming interface (API) that is
designed to help communicate the variables of interest to a user interface layer,
or possibly several different kinds of user interface layers at the same time,
without directly impacting how the models are constructed.

There are two different types of users for this API:

* **Model developers** should call :func:`export_variables` in their `build`
  methods fto tell the UI which variables are intended for the user.
* **User interface developers** should implement a module that creates a
  :class:`FlowsheetInterface` object and configures it with a flowsheet block,
  and appropriate methods bound to actions like `build` and `solve`.

The rest of this page will provide more detail for each type of user.

Model Developer Usage
---------------------

For a model developer, the primary interface is :func:`export_variables`. This function is applied to a
Pyomo block (also an IDAES one) to list the names (and, optionally some additional
information) of the variables that should be "exported" to the user interface.
It is expected that a standard set of exported variables will be performed by each
block (unit model, etc.) independently in the ``build`` method. For example, the
last line of the zero-order feed's `build` method is::

    export_variables(self, name="Feed Z0", desc="Zero-Order feed block",
                         variables=["flow_vol", "conc_mass_comp"])

At the time of that call, ``self`` is the feed block that was just built. It is
saying to export two variables called "flow_vol" and "conc_mass_comp".

.. autofunction:: export_variables

User Interface Developer Usage
------------------------------

There are two steps for creating a user interface to a flowsheet:

1. Define a function ``flowsheet_interface`` that creates the FlowsheetInterface.
   This object is not yet connected to an IDAES flowsheet
   block. The function should:

   a. Provide metadata (a name, description) for the flowsheet in the `FlowsheetInterface` constructor.

   b. Optionally list flowsheet-level variables (see :ref:`ui-define-variables`)

   c. Set the functions to call for "actions" -- by default, `build` and `solve` --
      that the UI can perform on the flowsheet.

   For example::

        def flowsheet_interface():
            fsi = FlowsheetInterface({"display_name": "METAB treatment train", "variables": []})
            fsi.set_action(WorkflowActions.build, build_flowsheet)
            fsi.set_action(WorkflowActions.solve, solve_flowsheet)
            return fsi

2. Define functions for the actions defined in Step 1. These fuunctions all have
   the following signature:

.. function:: action_function(block=None, ui=None, **kwargs)

    Perform an action.

    Args:
       block: Flowsheet block
       ui: FlowsheetInterface instance
       kwargs: Additional key

..

   For example::

        def build_flowsheet(ui=None, **kwargs):
            model = metab.build()  # 'metab' is the name of the flowsheet module
            # ..continue to build model..
            # last line of build should always be:
            ui.set_block(model)

        def solve_flowsheet(block=None, **kwargs):
            model = block
            metab.initialize_system(model)
            results = metab.solve(model)
            metab.assert_optimal_termination(results)


