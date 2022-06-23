Evaporation Pond (ZO)
=====================

Model Type
----------
This unit model is formulated as a single-input, double-output model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the cost_evaporation_pond method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Air temperature", "air_temperature", "K"
   "Daily solar radiation incident", "solar_radiation", "mJ/m**2"
   "Pond dike height", "dike_height", "ft"
   "Factor to adjust evaporation rate of pure water", "evaporation_rate_adj_factor", "None"
   "Evaporation rate calculation parameter A", "evap_rate_calc_a_parameter", "mm/d"
   "Evaporation rate calculation parameter B", "evap_rate_calc_b_parameter", "m**2/mJ"
   "Evaporation rate calculation parameter C", "evap_rate_calc_c_parameter", "m**2/mJ"
   "Adjusted area calculation parameter A", "adj_area_calc_a_parameter", "acre"
   "Adjusted area calculation parameter B", "adj_area_calc_b_parameter", "None"
   "Pond area needed based on evaporation rate", "area", "acre"
   "Adjusted pond area needed", "adj_area", "acre"
   "Calculated evaporation rate of pure water", "evaporation_rate_pure", "mm/d"
   "Pure water evaporation rate adjusted for salinity", "evaporation_rate_salt", "gal/acre/min"

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
