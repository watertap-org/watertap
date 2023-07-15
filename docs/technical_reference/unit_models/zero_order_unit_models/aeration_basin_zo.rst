Aeration Basin (ZO)
===================

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

.. index::
   pair: watertap.unit_models.zero_order.aeration_basin_zo;aeration_basin_zo

.. currentmodule:: watertap.unit_models.zero_order.aeration_basin_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.aeration_basin_zo
    :members:
    :noindex:
