CANDO-P (ZO)
============

Model Type
----------
This unit model is formulated as a **reactive single-inlet, double-outlet** model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.CANDOP_zo.CANDOPData.cost_CANDOP` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Electricity consumption of unit", "electricity", ":math:`kW`"
   "Electricity demand per kg N reacted", "electricity_intensity_N", ":math:`kWh/kg`"
   "Oxygen demand", "O2_demand", ":math:`kg/s`"
   "Oxygen consumed - nitrogen reacted ratio", "oxygen_nitrogen_ratio", ":math:`dimensionless`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for electricity consumption based on nitrogen consumed.", "electricity_consumption"
   "Constraint for oxygen consumption.", "oxygen_consumption"

.. index::
   pair: watertap.unit_models.zero_order.CANDOP_zo;CANDOP_zo

.. currentmodule:: watertap.unit_models.zero_order.CANDOP_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.CANDOP_zo
    :members:
    :noindex:
