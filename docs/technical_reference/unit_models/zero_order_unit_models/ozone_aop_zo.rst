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
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.ozone_aop_zo.OzoneAOPZOData.cost_ozonation_aop` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Ozone contact time", "contact_time", ":math:`min`"
   "CT value for ozone contactor", "concentration_time", ":math:`mg*min/l`"
   "Ozone mass transfer efficiency", "mass_transfer_efficiency", ":math:`dimensionless`"
   "Specific energy consumption for ozone generation", "specific_energy_coeff", ":math:`kWh/lb`"
   "Mass flow rate of ozone", "ozone_flow_mass", ":math:`lb/h`"
   "Ozone consumption", "ozone_consumption", ":math:`mg/l`"
   "Ozone generation power demand", "electricity", ":math:`kW`"
   "Oxidant dosage", "oxidant_dose", ":math:`mg/l`"
   "Mass flow rate of oxidant solution", "chemical_flow_mass", ":math:`kg/s`"
   "Ratio of ozone to total organic carbon", "ozone_toc_ratio", ":math:`dimensionless`"
   "Ratio of oxidant to ozone", "oxidant_ozone_ratio", ":math:`dimensionless`"

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
