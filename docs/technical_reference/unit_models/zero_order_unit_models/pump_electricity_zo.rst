Pump Electricity (ZO)
=====================

Model Type
----------
This unit model is formulated as a **pass-through** model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.pump_electricity_zo.PumpElectricityZOData.cost_pump_electricity` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Lift height for pump", "lift_height", ":math:`m`"
   "Efficiency of pump", "eta_pump", ":math:`dimensionless`"
   "Efficiency of motor", "eta_motor", ":math:`dimensionless`"
   "Electricity for low pressure pump", "electricity", ":math:`kW`"
   "Applied pressure", "applied_pressure", ":math:`bar`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for electricity consumption based on pump flowrate.", "electricity_consumption"
   "Constraint for pump applied pressure", "applied_pressure_constraint"

.. index::
   pair: watertap.unit_models.zero_order.pump_electricity_zo;pump_electricity_zo

.. currentmodule:: watertap.unit_models.zero_order.pump_electricity_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.pump_electricity_zo
    :members:
    :noindex:
