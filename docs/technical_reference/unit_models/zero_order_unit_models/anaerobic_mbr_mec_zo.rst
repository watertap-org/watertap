Integrated Anaerobic Membrane Bioreactor/Microbial Electrolysis Cell (ZO)
=========================================================================

Model Type
----------
This unit model is formulated as a **reactive single-inlet, double-outlet** model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **pump_electricity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.anaerobic_mbr_mec_zo.AnaerobicMBRMECZOData.cost_anaerobic_mbr_mec` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

.. index::
   pair: watertap.unit_models.zero_order.anaerobic_mbr_mec_zo;anaerobic_mbr_mec_zo

.. currentmodule:: watertap.unit_models.zero_order.anaerobic_mbr_mec_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.anaerobic_mbr_mec_zo
    :members:
    :noindex:
