Iron And Manganese Removal (ZO)
===============================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.iron_and_manganese_removal_zo.IronManganeseRemovalZOData.cost_iron_and_manganese_removal` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Ratio of air to water", "air_water_ratio", ":math:`dimensionless`"
   "Flow basis", "flow_basis", ":math:`m^3/h`"
   "Air flow rate", "air_flow_rate", ":math:`m^3/h`"
   "Constant in electricity intensity equation", "electricity_intensity_parameter", ":math:`h*hp/m^3`"
   "Dual media filter surface area", "filter_surf_area", ":math:`m^2`"
   "Number of dual media filter units", "num_filter_units", ":math:`dimensionless`"
   "Power consumption of iron and manganese removal", "electricity", ":math:`kW`"
   "Specific energy consumption with respect to feed flowrate", "electricity_intensity", ":math:`kWh/m^3`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Air flow rate constraint", "air_flow_rate_constraint"
   "Electricity intensity constraint", "electricity_intensity_constraint"
   "Power consumption constraint", "electricity_constraint"

.. index::
   pair: watertap.unit_models.zero_order.iron_and_manganese_removal_zo;iron_and_manganese_removal_zo

.. currentmodule:: watertap.unit_models.zero_order.iron_and_manganese_removal_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.iron_and_manganese_removal_zo
    :members:
    :noindex:
