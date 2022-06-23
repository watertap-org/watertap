Iron And Manganese Removal (ZO)
===============================

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
Costing is calculated using the cost_iron_and_manganese_removal method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

..csv-table::
   :header: "Description", "Variable Name"
   "Ratio of air to water", "air_water_ratio"
   "Flow basis", "flow_basis"
   "Air flow rate", "air_flow_rate"
   "Constant in electricity intensity equation", "electricity_intensity_parameter"
   "Dual media filter surface area", "filter_surf_area"
   "Number of dual media filter units", "num_filter_units"
   "Power consumption of iron and manganese removal", "electricity"
   "Specific energy consumption with respect to feed flowrate", "electricity_intensity"

.. index::
   pair: watertap.unit_models.zero_order.iron_and_manganese_removal_zo;iron_and_manganese_removal_zo

.. currentmodule:: watertap.unit_models.zero_order.iron_and_manganese_removal_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.iron_and_manganese_removal_zo
    :members:
    :noindex:
