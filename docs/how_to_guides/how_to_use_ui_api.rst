.. _howto_ui-api:

How-to guide for the UI API
===========================
.. py:currentmodule:: watertap.ui.fsapi

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

User interface developers can "discover" flowsheets for inclusion in the UI, and use the API to serialize and update from serialized data.

The rest of this page will provide more detail on how to do tasks for each type of user.

See also: :ref:`the UI flowsheet API reference section <ref_ui-fsapi>`.


Model Developers
----------------

.. image:: /_static/terminal-icon.png
    :height: 30px
    :align: left

Add an interface to your flowsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some Python module, define the function ``export_to_ui``, which will look
similar to this::

   def export_to_ui():
       return FlowsheetInterface(
           do_build=build_flowsheet,
           do_export=export_variables,
           do_solve=solve_flowsheet,
           name="My Flowsheet")

See :class:`FlowsheetInterface` for details on the arguments.

User Interface Developers
--------------------------

.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

.. _howto_api-find:

Find flowsheets
^^^^^^^^^^^^^^^^
Use the method :meth:`FlowsheetInterface.find` to get a mapping of module names to functions
that, when called, will create the flowsheet interface::

   results = fsapi.FlowsheetInterface.find("watertap")

Note that the returned interface will not create the flowsheet object and export the variables until the ``build`` method is invoked::

    first_module = list(results.keys())[0]
    interface = results[first_module]
    # at this point the name and description of the flowsheet are available
    interface.build()
    # at this point the flowsheet is built and all variables exported


.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

.. _howto_api-serialize:

Serialize flowsheet
^^^^^^^^^^^^^^^^^^^^
Use the ``dict()`` method to serialize the flowsheet::

    data = flowsheet.dict()

.. image:: /_static/menu-icon.png
    :height: 22px
    :align: left

.. _howto_api-load:

Load a serialized flowsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use the ``load()`` method to load values from a serialized flowsheet.
This will raise a ``MissingObjectError`` if any of the incoming values are not found in the target flowsheet::

   try:
       flowsheet.load(data)
   except fsapi.MissingObjectError as err:
       print(f"Error loading data: {err}")
       # print contents of list of missing items (key and variable name)
       for item in err.missing:
           print(f"Missing item: key={item.key}, name={item.name}")
