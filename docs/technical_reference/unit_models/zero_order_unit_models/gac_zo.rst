Granular Activated Carbon (ZO)
==============================

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
Costing is calculated using the cost_gac method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Empty bed contact time of unit", "empty_bed_contact_time", "hr"
   "Electricity consumption of unit", "electricity", "kW"
   "Parameter for calculating electricity based on empty bed contact time", "electricity_intensity_parameter", "kW/m**3"
   "Electricity intensity with respect to inlet flowrate of unit", "energy_electric_flow_vol_inlet", "kWh/m**3"
   "Replacement rate of activated carbon", "activated_carbon_replacement", "kg/m**3"
   "Demand for activated carbon", "activated_carbon_demand", "kg/hr"

Additional Constraints
----------------------

.. index::
   pair: watertap.unit_models.zero_order.gac_zo;gac_zo

.. currentmodule:: watertap.unit_models.zero_order.gac_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.gac_zo
    :members:
    :noindex:
