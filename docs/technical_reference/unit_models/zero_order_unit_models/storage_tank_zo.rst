Storage Tank (ZO)
=================

Model Type
----------
This unit model is formulated as a **pass-through** model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **constant_intensity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.storage_tank_zo.StorageTankZOData.cost_storage_tank` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Storage time needed", "storage_time", ":math:`h`"
   "Additional capacity needed for surge flow", "surge_capacity", ":math:`dimensionless`"
   "Storage tank volume", "tank_volume", ":math:`m^3`"

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
