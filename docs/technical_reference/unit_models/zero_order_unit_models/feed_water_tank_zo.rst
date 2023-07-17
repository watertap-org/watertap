Feed Water Tank (ZO)
====================

Model Type
----------
This unit model is formulated as a **pass-through** model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **constant_intensity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.core.zero_order_base.ZeroOrderBaseData.cost_power_law_flow` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

.. index::
   pair: watertap.unit_models.zero_order.feed_water_tank_zo;feed_water_tank_zo

.. currentmodule:: watertap.unit_models.zero_order.feed_water_tank_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.feed_water_tank_zo
    :members:
    :noindex:
