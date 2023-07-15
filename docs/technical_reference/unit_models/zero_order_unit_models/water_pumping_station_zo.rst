Water Pumping Station (ZO)
==========================

Model Type
----------
This unit model is formulated as a **pass-through** model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.core.zero_order_base.ZeroOrderBaseData.cost_power_law_flow` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Electricity consumption of unit", "electricity", ":math:`kW`"
   "Lift height for pump", "lift_height", ":math:`ft`"
   "Efficiency of pump", "eta_pump", ":math:`dimensionless`"
   "Efficiency of motor", "eta_motor", ":math:`dimensionless`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for electricity consumption based on pump flowrate.", "electricity_consumption"

.. index::
   pair: watertap.unit_models.zero_order.water_pumping_station_zo;water_pumping_station_zo

.. currentmodule:: watertap.unit_models.zero_order.water_pumping_station_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.water_pumping_station_zo
    :members:
    :noindex:
