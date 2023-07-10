Membrane Aerated Biofilm Reactor (ZO)
=====================================

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
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.mabr_zo.MABRZOData.cost_mabr` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Nitrogen removal rate per day", "nitrogen_removal_rate", ":math:`g/d/m^2`"
   "Sizing variable for effective reactor area", "reactor_area", ":math:`m^2`"
   "Air flow rate per area", "air_flow_rate", ":math:`m/h`"
   "Volumetric air flow rate", "air_flow_vol", ":math:`m^3/h`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for effective reactor area", "reactor_area_constraint"
   "Constraint for air flow", "air_flow_constraint"

.. index::
   pair: watertap.unit_models.zero_order.mabr_zo;mabr_zo

.. currentmodule:: watertap.unit_models.zero_order.mabr_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.mabr_zo
    :members:
    :noindex:
