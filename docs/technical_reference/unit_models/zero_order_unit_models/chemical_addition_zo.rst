Chemical Addition (ZO)
======================

Model Type
----------
This unit model is formulated as a pass-through model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the f(x) helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the cost_chemical_addition method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Dosing rate of chemical", "chemical_dosage", "mg/l"
   "Mass density of chemical solution", "solution_density", "kg/m**3"
   "Mass fraction of chemical in solution", "ratio_in_solution", "None"
   "Volumetric flow rate of chemical solution", "chemical_flow_vol", "m**3/s"
   "Electricity consumption of unit", "electricity", "kW"

Additional Constraints
----------------------

.. index::
   pair: watertap.unit_models.zero_order.chemical_addition_zo;chemical_addition_zo

.. currentmodule:: watertap.unit_models.zero_order.chemical_addition_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.chemical_addition_zo
    :members:
    :noindex:
