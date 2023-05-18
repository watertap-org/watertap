Gas Sparged Membrane (ZO)
=========================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **pump_electricity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
This unit does not include costing.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Mass flow of gas extracted per mass flow of influent", "gas_mass_influent_ratio", ":math:`dimensionless`"
   "Mass flow of hydrogen extracted", "flow_mass_gas_extraction", ":math:`kg/s`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Overall flow balance", "mass_balance"
   "Gas extraction equation", "mass_gas_extraction_equation"

.. index::
   pair: watertap.unit_models.zero_order.gas_sparged_membrane_zo;gas_sparged_membrane_zo

.. currentmodule:: watertap.unit_models.zero_order.gas_sparged_membrane_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.gas_sparged_membrane_zo
    :members:
    :noindex:
