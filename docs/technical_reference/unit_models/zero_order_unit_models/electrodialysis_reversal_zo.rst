Electrodialysis Reversal (ZO)
=============================

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
Costing is calculated using the cost_power_law_flow method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name"

   "Constant 1 in electricity intensity equation", "elec_coeff_1"
   "Constant 2 in electricity intensity equation", "elec_coeff_2"
   "Power consumption of brine concentrator", "electricity"
   "Specific energy consumption with respect to feed flowrate", "electricity_intensity"

.. index::
   pair: watertap.unit_models.zero_order.electrodialysis_reversal_zo;electrodialysis_reversal_zo

.. currentmodule:: watertap.unit_models.zero_order.electrodialysis_reversal_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.electrodialysis_reversal_zo
    :members:
    :noindex:
