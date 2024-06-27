Anaerobic Digestion Reactive (ZO)
=================================

Model Type
----------
This unit model is formulated as a **reactive single-inlet, double-outlet** model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

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

   "Ratio of m^3 biogas produced / kg TSS in influent", "biogas_tss_ratio", ":math:`m^3/kg`"
   "Biogas production", "biogas_production", ":math:`m^3/s`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for biogas production", "biogas_prod"

.. index::
   pair: watertap.unit_models.zero_order.anaerobic_digestion_reactive_zo;anaerobic_digestion_reactive_zo

.. currentmodule:: watertap.unit_models.zero_order.anaerobic_digestion_reactive_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.anaerobic_digestion_reactive_zo
    :members:
    :noindex:
