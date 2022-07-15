Ozone with Advanced Oxidation Processes (ZO)
============================================

Model Type
----------
This unit model is formulated as a **single-input, single-output** model form.
See documentation for :ref:`single-input, single-output Helper Methods<siso_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the **cost_ozonation_aop** method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Ozone contact time", "contact_time", "min"
   "CT value for ozone contactor", "concentration_time", "mg*min/l"
   "Ozone mass transfer efficiency", "mass_transfer_efficiency", "dimensionless"
   "Specific energy consumption for ozone generation", "specific_energy_coeff", "kWh/lb"
   "Mass flow rate of ozone", "ozone_flow_mass", "lb/hr"
   "Ozone consumption", "ozone_consumption", "mg/l"
   "Ozone generation power demand", "electricity", "kW"
   "Oxidant dosage", "oxidant_dose", "mg/l"
   "Mass flow rate of oxidant solution", "chemical_flow_mass", "kg/s"
   "Ratio of ozone to total organic carbon", "ozone_toc_ratio", "dimensionless"
   "Ratio of oxidant to ozone", "oxidant_ozone_ratio", "dimensionless"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Ozone consumption constraint", "ozone_consumption_constraint"
   "Ozone mass flow constraint", "ozone_flow_mass_constraint"
   "Ozone power constraint", "electricity_constraint"
   "Ozone/TOC ratio constraint", "ozone_toc_ratio_constraint"
   "Oxidant dose constraint", "oxidant_dose_constraint"
   "Oxidant mass flow constraint", "chemical_flow_mass_constraint"

.. index::
   pair: watertap.unit_models.zero_order.ozone_aop_zo;ozone_aop_zo

.. currentmodule:: watertap.unit_models.zero_order.ozone_aop_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.ozone_aop_zo
    :members:
    :noindex:
