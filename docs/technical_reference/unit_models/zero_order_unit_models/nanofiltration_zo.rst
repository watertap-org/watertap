Nanofiltration  (ZO)
====================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.nanofiltration_zo.NanofiltrationZOData.cost_nanofiltration` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Electricity consumption of unit", "electricity", ":math:`kW`"
   "Electricity intensity with respect to inlet flowrate of unit", "energy_electric_flow_vol_inlet", ":math:`kWh/m^3`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for electricity consumption based on feed flowrate.", "electricity_consumption"

.. index::
   pair: watertap.unit_models.zero_order.nanofiltration_zo;nanofiltration_zo

.. currentmodule:: watertap.unit_models.zero_order.nanofiltration_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.nanofiltration_zo
    :members:
    :noindex:
