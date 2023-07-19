Filter Press (ZO)
=================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.filter_press_zo.FilterPressZOData.cost_filter_press` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Hours per day of filter press operation", "hours_per_day_operation", ":math:`h/d`"
   "Filter press cycle time", "cycle_time", ":math:`h`"
   "Parameter A for electricity calculation", "electricity_a_parameter", ":math:`kWh/a/ft^3`"
   "Parameter B for electricity calculation", "electricity_b_parameter", ":math:`dimensionless`"
   "Filter press capacity", "filter_press_capacity", ":math:`ft^3`"
   "Filter press power", "electricity", ":math:`kW`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Filter press capacity constraint", "fp_capacity"
   "Filter press electricity constraint", "fp_electricity"

.. index::
   pair: watertap.unit_models.zero_order.filter_press_zo;filter_press_zo

.. currentmodule:: watertap.unit_models.zero_order.filter_press_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.filter_press_zo
    :members:
    :noindex:
