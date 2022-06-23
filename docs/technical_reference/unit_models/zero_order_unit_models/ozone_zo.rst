Ozone  (ZO)
===========

Model Type
----------
This unit model is formulated as a single-input, single-output model form.
See documentation for :ref:`single-input, single-output Helper Methods<siso_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the f(x) helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the cost_ozonation method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Ozone contact time", "contact_time", "min"
   "CT value for ozone contactor", "concentration_time", "mg*min/l"
   "Ozone mass transfer efficiency", "mass_transfer_efficiency", "None"
   "Specific energy consumption for ozone generation", "specific_energy_coeff", "kWh/lb"
   "Mass flow rate of ozone", "ozone_flow_mass", "lb/hr"
   "Ozone consumption", "ozone_consumption", "mg/l"
   "Ozone generation power demand", "electricity", "kW"

Additional Constraints
----------------------

.. index::
   pair: watertap.unit_models.zero_order.ozone_zo;ozone_zo

.. currentmodule:: watertap.unit_models.zero_order.ozone_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.ozone_zo
    :members:
    :noindex:
