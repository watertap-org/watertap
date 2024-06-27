Microbial Battery (ZO)
======================

Model Type
----------
This unit model is formulated as a **reactive single-inlet, double-outlet** model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **constant_intensity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.microbial_battery_zo.MicrobialBatteryData.cost_microbial_battery` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Hydraulic retention time of water treatment unit", "HRT", ":math:`h`"
   "Volume of water treatment unit", "reactor_volume", ":math:`m^3`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for reactor volume.", "reactor_volume_rule"

.. index::
   pair: watertap.unit_models.zero_order.microbial_battery_zo;microbial_battery_zo

.. currentmodule:: watertap.unit_models.zero_order.microbial_battery_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.microbial_battery_zo
    :members:
    :noindex:
