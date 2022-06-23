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
   :header: "Description", "Variable Name"

   "Empty bed contact time of unit", "empty_bed_contact_time"
   "Electricity consumption of unit", "electricity"
   "Parameter for calculating electricity based on empty bed contact time", "electricity_intensity_parameter"
   "Electricity intensity with respect to inlet flowrate of unit", "energy_electric_flow_vol_inlet"
   "Replacement rate of activated carbon", "activated_carbon_replacement"
   "Demand for activated carbon", "activated_carbon_demand"

.. index::
   pair: watertap.unit_models.zero_order.gac_zo;gac_zo

.. currentmodule:: watertap.unit_models.zero_order.gac_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.gac_zo
    :members:
    :noindex:
