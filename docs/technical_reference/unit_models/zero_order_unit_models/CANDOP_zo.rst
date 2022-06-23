CANDO-P (ZO)
============

Model Type
----------
This unit model is formulated as a reactive single-inlet, double-outlet model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the f(x) helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the cost_CANDOP method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name"

   "Electricity consumption of unit", "electricity"
   "Electricity demand per kg N reacted", "electricity_intensity_N"
   "Oxygen demand", "O2_demand"
   "Oxygen consumed - nitrogen reacted ratio", "oxygen_nitrogen_ratio"

.. index::
   pair: watertap.unit_models.zero_order.CANDOP_zo;CANDOP_zo

.. currentmodule:: watertap.unit_models.zero_order.CANDOP_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.CANDOP_zo
    :members:
    :noindex:
