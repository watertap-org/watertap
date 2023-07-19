Volatile Fatty Acid (VFA) Recovery Unit (ZO)
============================================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **constant_intensity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.vfa_recovery_zo.VFARecoveryZOData.cost_vfa_recovery` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Thermal energy required per mass VFA", "heat_required_per_vfa_mass", ":math:`kJ/kg`"
   "Thermal energy required", "heat_consumption", ":math:`kJ/s`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for heat consumption", "eq_heat_consumption"

.. index::
   pair: watertap.unit_models.zero_order.vfa_recovery_zo;vfa_recovery_zo

.. currentmodule:: watertap.unit_models.zero_order.vfa_recovery_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.vfa_recovery_zo
    :members:
    :noindex:
