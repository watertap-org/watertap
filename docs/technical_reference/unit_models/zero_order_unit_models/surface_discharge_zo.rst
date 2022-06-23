Surface Discharge (ZO)
======================

Model Type
----------
This unit model is formulated as a pass-through model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the cost_surface_discharge method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Piping distance", "pipe_distance", "mi"
   "Pipe diameter", "pipe_diameter", "in"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"


.. index::
   pair: watertap.unit_models.zero_order.surface_discharge_zo;surface_discharge_zo

.. currentmodule:: watertap.unit_models.zero_order.surface_discharge_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.surface_discharge_zo
    :members:
    :noindex:
