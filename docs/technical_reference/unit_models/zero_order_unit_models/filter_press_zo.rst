Filter Press (ZO)
=================

Model Type
----------
This unit model is formulated as a single-input, double-output model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the f(x) helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the cost_filter_press method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name"

   "Hours per day of filter press operation", "hours_per_day_operation"
   "Filter press cycle time", "cycle_time"
   "Parameter A for electricity calculation", "electricity_a_parameter"
   "Parameter B for electricity calculation", "electricity_b_parameter"
   "Filter press capacity", "filter_press_capacity"
   "Filter press power", "electricity"

.. index::
   pair: watertap.unit_models.zero_order.filter_press_zo;filter_press_zo

.. currentmodule:: watertap.unit_models.zero_order.filter_press_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.filter_press_zo
    :members:
    :noindex:
