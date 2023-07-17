Recirculating Dynamic Membrane Bioreactor (ZO)
==============================================

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
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.dmbr_zo.DMBRZOOData.cost_dmbr` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

.. index::
   pair: watertap.unit_models.zero_order.dmbr_zo;dmbr_zo

.. currentmodule:: watertap.unit_models.zero_order.dmbr_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.dmbr_zo
    :members:
    :noindex:
