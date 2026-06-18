Crystallizer (ZO)
=================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.crystallizer_zo.CrystallizerZOData.cost_crystallizer` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Constant 1 in electricity intensity equation", "elec_coeff_1", ":math:`kWh/m^3`"
   "Constant 2 in electricity intensity equation", "elec_coeff_2", ":math:`kWh/m^3/(mg/L)`"
   "Constant 3 in electricity intensity equation", "elec_coeff_3", ":math:`kWh/m^3`"
   "Constant 4 in electricity intensity equation", "elec_coeff_4", ":math:`kWh/m^3/(m^3/h)`"
   "Power consumption of crystallizer", "electricity", ":math:`kW`"
   "Specific energy consumption with respect to feed flowrate", "electricity_intensity", ":math:`kWh/m^3`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Electricity intensity constraint", "electricity_intensity_constraint"
   "Power consumption constraint", "electricity_constraint"

.. index::
   pair: watertap.unit_models.zero_order.crystallizer_zo;crystallizer_zo

.. currentmodule:: watertap.unit_models.zero_order.crystallizer_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.crystallizer_zo
    :members:
    :noindex:
