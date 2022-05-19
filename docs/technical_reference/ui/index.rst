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

* :ref:`Model developers <ui-model-dev>` should call :func:`export_variables` in their `build`
  methods to tell the UI which variables are intended for the user.
* :ref:`User interface developers <ui-interface-dev>` should implement a module that creates a
  :class:`FlowsheetInterface` object and configures it with a flowsheet block,
  and appropriate methods bound to actions like `build` and `solve`.

The rest of this page will provide more detail for each type of user.

.. _ui-model-dev:

Model Developer Usage
---------------------

For a model developer, the primary interface is the function :func:`export_variables`.
This function is applied to a Pyomo block (also an IDAES one) to list the names (and, optionally some additional information) of the variables that should be "exported" to the user interface.
It is expected that a standard set of exported variables will be performed by each block (unit model, etc.) independently in the ``build`` method.
For example, the last line of the zero-order feed's `build` method is::

    export_variables(self, name="Feed Z0", desc="Zero-Order feed block",
                         variables=["flow_vol", "conc_mass_comp"])

At the time of that call, ``self`` is the feed block that was just built.
It is saying to export two variables called "flow_vol" and "conc_mass_comp".

Testing the model interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^
To make sure you are exporting the variables you meant to, you can retrieve them as a list.
This requires two calls:

1. Get the interface that is attached to the block by ``export_variables`` by calling :func:`get_block_interface`.
2. Iterate through variables exported by that block by calling :meth:`BlockInterface.get_exported_variables`.
   Each variable that is returned will have the same structure as an entry in the "variables" ConfigList of the :attr:`BlockInterface.CONFIG`.
   For example, if the model block was called ``model.fs.thing``, you could do the following::

        def test_variables_in_thing():
            model = somehow_create_the_model()
            expected_thing_names = ["foo", "bar"]
            for thing_var in get_block_interface(model.fs.thing).get_exported_variables():
                 assert thing_var["name"] in expected_thing_names




.. _ui-interface-dev:

User Interface Developer Usage
------------------------------
The user interface developer has to do two primary tasks:

* :ref:`Create interfaces <ui-create-interface>` to specific flowsheets
* :ref:`Find and use those interfaces <ui-finduse-interface>` in the logic of the UI backend

The following two sections (linked above) show how to perform those tasks.

.. _ui-create-interface:

Creating an interface to a flowsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two steps for **creating a user interface** to a flowsheet:

1. Define a function ``flowsheet_interface()`` that creates the :class:`FlowsheetInterface` object.
   This function should have the following signature:

.. function:: flowsheet_interface() -> FlowsheetInterface

   In other words, it takes no arguments and returns a :class:`FlowsheetInterface` object.
   This object is not yet connected to an IDAES flowsheet block.
   The function should (a) define the metadata and variables for the flowsheet using the ConfigDict documented in :attr:`BlockInterface.config`, then
   (b) set the functions to call for the "actions" -- by default, `build` and `solve` -- that the UI can perform on the flowsheet.

   For example::

        def flowsheet_interface():
            fsi = FlowsheetInterface({
              "name": "metab",
              "display_name": "METAB treatment train",
              "variables": [
                  {"name": "var1", "display_name": "Variable Numero Uno",
                   "description": "The first of the variables", "units": "m**3"},
                  {"name": "var2", "display_name": "Variable Numero Dos",
                   "description": "The second of the variables", "units": "m**4"}
              ]
            })
            fsi.set_action(WorkflowActions.build, build_flowsheet)
            fsi.set_action(WorkflowActions.solve, solve_flowsheet)
            return fsi

Note that you only need to add variables that are not already exported by the model, and that there are pretty reasonable defaults for things like the name, display_name (same as name), and description. So in most cases this will be a very simple call; the extended version was shown for didactic purposes.

2. Define functions for the actions defined in Step 1. These functions all have the following signature:

.. function:: action_function([block=None, ui=None], **kwargs)

    Perform an action.

    :param Block block: Flowsheet block
    :param FlowsheetInterface ui: FlowsheetInterface instance
    :param dict kwargs: Additional key/value pairs specific to this action

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


If you wish to define your own actions, use the :meth:`~FlowsheetInterface.add_action_type` method of the object that was created by ``flowsheet_interface()``.

.. _ui-define-variables:

Defining variables in more detail
+++++++++++++++++++++++++++++++++

.. todo: two ways to do it (1) provide more infor to export_variables, (2) create the BlockInterface yourself.

.. _ui-finduse-interface:

Find and use flowsheet interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have created a flowsheet interface, as described in :ref:`ui-create-interface`, you need to use it in the UI backend.
