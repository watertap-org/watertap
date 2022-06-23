Brine Concentrator (ZO)
=======================

Model Type
----------
This unit model is formulated as a single-input, double-output model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the f(x) helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the cost_brine_concentrator method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Constant 1 in electricity intensity equation", "elec_coeff_1", "kWh/m**3"
   "Constant 2 in electricity intensity equation", "elec_coeff_2", "kWh*l/m**3/mg"
   "Constant 3 in electricity intensity equation", "elec_coeff_3", "kWh/m**3"
   "Constant 4 in electricity intensity equation", "elec_coeff_4", "hr*kWh/m**6"
   "Power consumption of brine concentrator", "electricity", "kW"
   "Specific energy consumption with respect to feed flowrate", "electricity_intensity", "kWh/m**3"

Additional Constraints
----------------------

.. index::
   pair: watertap.unit_models.zero_order.brine_concentrator_zo;brine_concentrator_zo

.. currentmodule:: watertap.unit_models.zero_order.brine_concentrator_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.brine_concentrator_zo
    :members:
    :noindex:
