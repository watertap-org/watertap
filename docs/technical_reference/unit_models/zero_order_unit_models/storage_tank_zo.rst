Storage Tank (ZO)
=================

Model Type
----------
This unit model is formulated as a pass-through model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the cost_storage_tank method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Storage time needed", "storage_time", "hr"
   "Additional capacity needed for surge flow", "surge_capacity", "None"
   "Storage tank volume", "tank_volume", "m**3"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Tank volume constraint", "tank_volume_constraint"

.. index::
   pair: watertap.unit_models.zero_order.storage_tank_zo;storage_tank_zo

.. currentmodule:: watertap.unit_models.zero_order.storage_tank_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.storage_tank_zo
    :members:
    :noindex:
