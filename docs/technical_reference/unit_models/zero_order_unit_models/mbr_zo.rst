Membrane Bioreactor (ZO)
========================

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
Costing is calculated using the :py:meth:`~watertap.core.zero_order_base.ZeroOrderBaseData.cost_power_law_flow` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Constant 1 in electricity intensity equation", "elec_coeff_1", ":math:`kWh/m^3`"
   "Constant 2 in electricity intensity equation", "elec_coeff_2", ":math:`dimensionless`"
   "Specific energy consumption with respect to feed flowrate", "electricity_intensity", ":math:`kWh/m^3`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Electricity intensity constraint", "electricity_intensity_constraint"
   "Power consumption constraint", "electricity_constraint"

.. index::
   pair: watertap.unit_models.zero_order.mbr_zo;mbr_zo

.. currentmodule:: watertap.unit_models.zero_order.mbr_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.mbr_zo
    :members:
    :noindex:
