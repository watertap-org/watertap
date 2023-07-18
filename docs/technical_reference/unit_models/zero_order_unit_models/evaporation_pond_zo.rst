Evaporation Pond (ZO)
=====================

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
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.evaporation_pond_zo.EvaporationPondZOData.cost_evaporation_pond` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Air temperature", "air_temperature", ":math:`K`"
   "Daily solar radiation incident", "solar_radiation", ":math:`mJ/m^2`"
   "Pond dike height", "dike_height", ":math:`ft`"
   "Factor to adjust evaporation rate of pure water", "evaporation_rate_adj_factor", ":math:`dimensionless`"
   "Evaporation rate calculation parameter A", "evap_rate_calc_a_parameter", ":math:`mm/d`"
   "Evaporation rate calculation parameter B", "evap_rate_calc_b_parameter", ":math:`m^2/mJ`"
   "Evaporation rate calculation parameter C", "evap_rate_calc_c_parameter", ":math:`m^2/mJ`"
   "Adjusted area calculation parameter A", "adj_area_calc_a_parameter", ":math:`acre`"
   "Adjusted area calculation parameter B", "adj_area_calc_b_parameter", ":math:`dimensionless`"
   "Pond area needed based on evaporation rate", "area", ":math:`acre`"
   "Adjusted pond area needed", "adj_area", ":math:`acre`"
   "Calculated evaporation rate of pure water", "evaporation_rate_pure", ":math:`mm/d`"
   "Pure water evaporation rate adjusted for salinity", "evaporation_rate_salt", ":math:`gal/acre/min`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Evaporation rate of pure water constraint", "evap_rate_pure_constraint"
   "Adjusted evaporation rate for salinity constraint", "evap_rate_salt_constraint"
   "Base area constraint", "area_constraint"
   "Adjusted area constraint", "area_adj_constraint"

.. index::
   pair: watertap.unit_models.zero_order.evaporation_pond_zo;evaporation_pond_zo

.. currentmodule:: watertap.unit_models.zero_order.evaporation_pond_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.evaporation_pond_zo
    :members:
    :noindex:
