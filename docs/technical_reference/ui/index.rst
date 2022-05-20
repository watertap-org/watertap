WaterTAP User Interface API
===========================

.. include:: <isonum.txt>

.. py:currentmodule:: watertap.ui.api

.. contents:: Contents
    :depth: 1
    :local:

Introduction
------------

This page describes an application programming interface (API) that is designed to help communicate information about a flowsheet and its variables to a user interface layer (or possibly several different kinds of user interface layers at the same time).
The API also provides the ability, intended for the UI developer, to update the variable values and run "actions" such as building and solving the flowsheet.

There are two distinct intended users for this API:

.. image:: /_static/terminal-icon.png
    :height: 20px
    :align: left

:ref:`Model developers <ui-model-dev>` can select which variables to "export" to the UI for each component of the model, and provide extra metadata (display name, description) for them.

.. image:: /_static/menu-icon.png
    :height: 15px
    :align: left

:ref:`User interface developers <ui-interface-dev>` can select a flowsheet to wrap, and implement the interface to the build, solve, and possibly other actions that the UI could take with the flowsheet.

The rest of this page will provide more detail for each type of user.

.. image:: /_static/terminal-icon.png
    :height: 60px
    :align: left

.. _ui-model-dev:

Model Developer Usage
---------------------

|


For a model developer, the primary interface is the function :func:`export_variables`.
This function is applied to a Pyomo block (also an IDAES one) to list the names (and, optionally some additional information) of the variables that should be "exported" to the user interface.
It is expected that a standard set of exported variables will be performed by each block (unit model, etc.) independently in the ``build`` method. To get all variables exported by the flowsheet, all the exports of sub-blocks will be gathered into a hierarchy.
For example, the last line of the zero-order feed's `build` method is::

    export_variables(self, name="Feed Z0", desc="Zero-Order feed block",
                         variables=["flow_vol", "conc_mass_comp"])

At the time of that call, ``self`` is the feed block that was just built.
It is saying to export two variables called "flow_vol" and "conc_mass_comp".

There are other things about variables that can be specified, in which case instead of the simple variable name you need to use a dictionary with fields given by the "variables" section of the :class:`BlockInterface` configuration options.

For example, if you wanted to mark "flow_vol" as read-only (loading it into the model will not change it), then you could do::

    export_variables(self, name="Feed Z0", desc="Zero-Order feed block",
                         variables=[{"name":"flow_vol", "readonly": True},
                                    "conc_mass_comp"])



Lower-level interface
^^^^^^^^^^^^^^^^^^^^^
Under the covers, the ``export_variables`` function creates an instance of :class:`BlockInterface` bound to the given IDAES model component, and configures it metadata for itself and its variables.
Therefore, you can directly call the ``BlockInterface`` constructor with the block object and a dictionary for configuration options.

For example, for the `export_variables` call shown above, the equivalent lower-level calls would look like this::

    # as above, assume the model is in "self"
    BlockInterface(self, {
        name: "Feed Z0",
        description: "Zero-Order feed block",
        "variables": [
            {"name": "flow_vol"},
            {"name": "conc_mass_comp"}
        ]
    })



.. image:: /_static/menu-icon.png
    :height: 65px
    :align: left
.. _ui-interface-dev:

User Interface Developer Usage
------------------------------

|

The user interface developer will use the Python API for two tasks:

* :ref:`Create interfaces <ui-create-interface>` to specific flowsheets
* :ref:`Find and use those interfaces <ui-finduse-interface>` in the web server (or other UI backend)

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


If you wish to define your own actions for a given flowsheet, use the :meth:`~FlowsheetInterface.add_action_type` method.

Action dependencies
+++++++++++++++++++
You are able to specify dependencies of an action on other actions, which means that this action will automatically run the other actions (if they have not been already run).

The execution of dependencies are tracked so that they have the behavior one would expect, and are not re-run unless necessary. For example, for the built-in actions, where
**get-results** depends on **solve** and
**solve** depends on **build**, the following sequence of calls would run actions as shown:

1. **solve** |rarr| run `build`, `solve`
2. **solve** |rarr| run `solve` (`build` is not re-run)
3. **get-results** |rarr| run `get-results`
4. **build** |rarr| run `build`
5. **get-results** |rarr| run `solve`, `get-results` (running `build` in (4) means `solve` needs to be re-run)


.. _ui-finduse-interface:

Find and use flowsheet interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. contents:: Contents
    :local:
    :depth: 2

Once you have created a flowsheet interface, as described in :ref:`ui-create-interface`, you need to use it in the UI backend.


.. TODO Implement on backend, then return and document here.

Finding flowsheets
++++++++++++++++++
**TBD**

Interacting with flowsheets
+++++++++++++++++++++++++++
**TBD**

Fetch and update flowsheet values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**TBD**

Get the flowsheet diagram
~~~~~~~~~~~~~~~~~~~~~~~~~
**TBD**

Run flowsheet actions
~~~~~~~~~~~~~~~~~~~~~
**TBD**

Save flowsheet status in a file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**TBD**

Flowsheet information format
++++++++++++++++++++++++++++
The information about the flowsheet and all its subblocks and variables is encoded in JSON when it is transferred between the backend and the UI, or saved to a file. The form of this information is shown in the following JSON schema:

.. image:: /_static/search-icon.png
    :height: 65px
    :align: left

Class and function reference
----------------------------

|

High-level API
^^^^^^^^^^^^^^

.. autofunction:: export_variables

.. autoclass:: FlowsheetInterface
    :members:
    :special-members: __init__, __eq__

Lower-level API
^^^^^^^^^^^^^^^

.. autofunction:: get_block_interface

.. autofunction:: set_block_interface

.. autoclass:: BlockInterface
    :members:
    :special-members: __init__

.. autoclass:: WorkflowActions
    :members:
    :undoc-members:
