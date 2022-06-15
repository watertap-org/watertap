.. _howto_ui-api:

How-to guide for the UI API
===========================
.. py:currentmodule:: watertap.ui.api

.. contents:: Contents
    :depth: 2
    :local:

Overview
--------

There are two distinct intended users for this API:

.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

Model developers  can select which variables to "export" to the UI for each component of the model, and provide extra metadata (display name, description) for them.
For flowsheets, they also should specify how to build and solve the flowsheets.

.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

User interface developers can "discover" flowsheets for inclusion in the UI, and use the API to provide JSON versions of their contents (and results).

The rest of this page will provide more detail on how to do tasks for each type of user.

See also: :ref:`the UI API reference section <ref_ui-api>`.


Model Developers
----------------

.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

Exporting variables with :func:`export_variables`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each IDAES component and/or Pyomo block the developer should list the names of the variables that should be "exported" to the user interface. 
Optionally, users may provide additional information as well (see :func:`export_variables` for details).
This should be done in the component's ``build`` method.
Each component/block need only worry about its own variables.
Teh API will take care of gathering all the components for a flowsheet together.
For example, the last line of the zero-order feed's `build` method is::

    export_variables(self, name="Feed Z0", desc="Zero-Order feed block",
                         variables=["flow_vol", "conc_mass_comp"])

At the time of that call, ``self`` is the feed block that was just built.
It is saying to export two variables called "flow_vol" and "conc_mass_comp".

There are other things about variables that can be specified.
For details see :func:`export_variables`.

For example, if you wanted to mark "flow_vol" as read-only (loading it into the model will not change it), then you could do::

    export_variables(self, name="Feed Z0", desc="Zero-Order feed block",
                         variables={"flow_vol": {"readonly": True},
                                    "conc_mass_comp": {}})


.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

.. _howto_api-ui_create-interface:

Creating an interface to a flowsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three steps for **creating a user interface** to a flowsheet:

1. :ref:`Define a function called flowsheet_interface() <howto_api-ui_define-interface>` that creates the :class:`FlowsheetInterface` object.
2. Define flowsheet actions. See :ref:`howto_api-ui_create-action`.
   All interfaces should define the following actions:

    build
        Build a flowsheet. At the end this will set the built flowsheet into the provided 'ui' parameter with ``ui.set_block(<flowsheet-that-was-built>)``.

    solve
        Solve a flowsheet

3. :ref:`Set the action functions into the flowsheet interface <howto_api-set-action>` before returning it to the user.


.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

.. _howto_api-ui_define-interface:

Define the ``flowsheet_interface`` function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The flowsheet interface is created by a function with this signature:

.. function:: flowsheet_interface() -> FlowsheetInterface

   In other words, it takes no arguments and returns a :class:`FlowsheetInterface` object.
   This object is not yet connected to an IDAES flowsheet block.
   The function should return a flowsheet interface, see :ref:`howto_api-ui_create-interface`.

Note that you only need to add variables that are not already exported by the model, and that there are reasonable defaults for things like the name, display_name (same as name), and description. So in most cases this will be a very simple call; the extended version was shown for didactic purposes.

.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

.. _howto_api-ui_create-action:

Create an "action" function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The action functions should operate on optional keywords for the flowsheet block and the FlowsheetInterface instance.
It will call the WaterTAP code to perform the appropriate actions.
The function should have the following signature:

.. function:: action_function([block=None, ui=None], **kwargs)

    Perform an action.

    :param Block block: Flowsheet block
    :param FlowsheetInterface ui: FlowsheetInterface instance
    :param dict kwargs: Additional key/value pairs specific to this action

For example::

    def build_flowsheet(ui=None, **kwargs):
        model = my_model.build()
        my_model.set_operating_conditions(model)
        my_model.assert_degrees_of_freedom(model, 0)
        my_model.assert_units_consistent(model)
        my_model.add_costing(model)
        model.fs.costing.initialize()
        # Export some additional costing variables
        export_variables(
            model.fs.costing,
            name="My model costing",
            desc="Costing block for METAB model",
            category="costing",
            variables=[
                "utilization_factor",
                "TIC",
                "maintenance_costs_percent_FCI",
            ],
        )
        my_model.adjust_default_parameters(model)
        my_model.assert_degrees_of_freedom(model, 0)

        # ** IMPORTANT **
        # Set this flowsheet as the top-level block for the interface
        ui.set_block(model.fs)



.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

.. _howto_api-set-action:

Set actions for a flowsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once you have :ref:`created a flowsheet action <howto_api-ui_create-action>`, you can set it in the flowsheet by calling :meth:`FlowsheetInterface.set_action()`. For example::

    from watertap.ui.api import FlowsheetInterface, WorkflowActions
    def flowsheet_interface():
        fsi = FlowsheetInterface({"display_name": "My treatment train",
                                  "description": "Treatment train to show off my model"})
        fsi.set_action(WorkflowActions.build, build_flowsheet)
        fsi.set_action(WorkflowActions.solve, solve_flowsheet)
        return fsi


User Interface Developers
--------------------------

.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

.. _howto_api-finduse-interface:

Find flowsheets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use the function :func:`find_flowsheet_interfaces` to get a list of :class:`FlowsheetInterface` objects representing the available interfaces.
By default this function will find all interfaces in the ``watertap`` Python package in which it is situated::

    from watertap.ui.api import find_flowsheet_interfaces
    for fsi in find_flowsheet_interfaces():
        print(f"Got flowsheet: {fsi.display_name}")

You can add a configuration file or dict to specify alternate or additional places to look.
Just remember when doing so to provide the default watertap location if you want it included.
For example::

    from watertap.ui.api import find_flowsheet_interfaces
    interface_list = find_flowsheet_interfaces(config={
        "packages": ["watertap", "my_other_package"]})
    for fsi in interface_list:
        print(f"Got flowsheet: {fsi.display_name}")


.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

.. _howto_api-fetch-update:

Fetch and update flowsheet values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The values for all variables exported by the flowsheet are available via the ``dict()`` method.
The format of the returned value is documented in the :mod:`~watertap.ui.api` module header.
Note: if you want to write these values as JSON to a stream, use :meth:`FlowsheetInterface.save`.
For example::

    from watertap.ui.api import find_flowsheet_interfaces
    for fsi in find_flowsheet_interfaces():
        print(f"Got flowsheet: {fsi.display_name}")
        print(f"Flowsheet contents: {fsi.dict()}")

.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

.. _howto_api-run-actions:

Run flowsheet actions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once flowsheet actions are created, you invoke them with the :meth:`FlowsheetInterface.run_action` method. The name of the action, if it is a standard one, will be an attribute in the :class:`WorkflowActions` class. For example::

    from watertap.ui.api import FlowsheetInterface

    def run_build(fsi: FlowsheetInterface):
        fsi.run_action(WorkflowActions.build)

.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

.. _howto_api-save-flowsheet:

Save flowsheet values in a file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The values of variables that were :func:`exported <export_variables>` for a given flowsheet can be saved to a file.
For example::

    from watertap.ui.api import FlowsheetInterface
    fsi = FlowsheetInterface()
    fsi.set_block(my_flowsheet)
    fsi.save("my_flowsheet_saved.json")


Note: The method is simply a wrapper that calls ``dict()`` and feeds the result to a JSON serializer.

.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

.. _howto_api-load-flowsheet:

Load flowsheet values from a file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The values of variables that were :ref:`saved <howto_api-save-flowsheet>` for a given flowsheet can be loaded back into the model.
This operation changes the matching values in the actual underlying model to the values that are stored in the file.
In other words, unlike the :meth:`FlowsheetInterface.save` method, this method changes what is stored in memory.
Invoking load is straightforward::

    from watertap.ui.api import FlowsheetInterface
    # Create and save
    fsi = FlowsheetInterface()
    fsi.set_block(my_flowsheet)
    fsi.save("my_flowsheet_saved.json")
    # Load back
    fsi.load("my_flowsheet_saved.json")

To handle situations where the model changes over time, and the saved data and model do not match exactly, there are two methods you can use after you are done loading:

    :meth:`FlowsheetInterface.get_var_missing`
        Returns variables that were in the loaded data, but not in the underlying model.

    :meth:`FlowsheetInterface.get_var_extra`
        Returns variables that are exported by the model, but were not in the loaded data.

In both cases, the method returns a mapping of with the full name of the block (e.g., `flowsheet.component.subcomponent`) as the key and a list of variable names as the value.