Coagulation Flocculation (ZO)
=============================

Model Type
----------
This unit model is formulated as a **pass-through** model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.coag_and_floc_zo.CoagulationFlocculationZOData.cost_coag_and_floc` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Dose of alum", "alum_dose", ":math:`mg/l`"
   "Mass ratio of alum in feed solution", "alum_ratio_in_solution", ":math:`dimensionless`"
   "Mass density of alum feed solution", "alum_solution_density", ":math:`kg/m^3`"
   "Mass flow rate of alum", "alum_flow_mass", ":math:`kg/s`"
   "Volumetric flow rate of alum solution", "alum_flow_vol", ":math:`m^3/s`"
   "Dose of polymer", "polymer_dose", ":math:`mg/l`"
   "Mass ratio of polymer in feed solution", "polymer_ratio_in_solution", ":math:`dimensionless`"
   "Mass density of polymer feed solution", "polymer_solution_density", ":math:`kg/m^3`"
   "Mass flow rate of polymer", "polymer_flow_mass", ":math:`kg/s`"
   "Volumetric flow rate of polymer solution", "polymer_flow_vol", ":math:`m^3/s`"
   "Pump head for chemical injection pumps", "chem_pump_head", ":math:`ft`"
   "Pump efficiency for chemical injection pumps", "eta_chem_pump", ":math:`dimensionless`"
   "Rapid Mix Retention Time", "rapid_mix_retention_time", ":math:`s`"
   "Floc Retention Time", "floc_retention_time", ":math:`min`"
   "Rapid Mix Basin Volume", "rapid_mix_basin_vol", ":math:`m^3`"
   "Floc Basin Volume", "floc_basin_vol", ":math:`m^3`"
   "Number of Rapid Mixers", "num_rapid_mixers", ":math:`dimensionless`"
   "Number of Floc Mixers", "num_floc_mixers", ":math:`dimensionless`"
   "Number of Rapid Mix Processes", "num_rapid_mix_processes", ":math:`dimensionless`"
   "Number of Floc Processes", "num_floc_processes", ":math:`dimensionless`"
   "Number of Coagulation Processes", "num_coag_processes", ":math:`dimensionless`"
   "Number of Floc Injection Processes", "num_floc_injection_processes", ":math:`dimensionless`"
   "Rapid Mix Velocity Gradient", "velocity_gradient_rapid_mix", ":math:`1/s`"
   "Floc Velocity Gradient", "velocity_gradient_floc", ":math:`1/s`"
   "Alum Addition Power Consumption", "power_alum_addition", ":math:`kW`"
   "Polymer Addition Power Consumption", "power_polymer_addition", ":math:`kW`"
   "Rapid Mix Power Consumption", "power_rapid_mix", ":math:`kW`"
   "Floc Power Consumption", "power_floc", ":math:`kW`"
   "Total Power Consumption", "electricity", ":math:`kW`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "None", "rapid_mix_basin_vol_constraint"
   "None", "floc_basin_vol_constraint"
   "None", "alum_flow_mass_constraint"
   "None", "alum_flow_vol_constraint"
   "None", "polymer_flow_mass_constraint"
   "None", "polymer_flow_vol_constraint"
   "None", "rapid_mix_power_constraint"
   "None", "floc_power_constraint"
   "None", "alum_addition_power_constraint"
   "None", "polymer_addition_power_constraint"
   "None", "electricity_constraint"

.. index::
   pair: watertap.unit_models.zero_order.coag_and_floc_zo;coag_and_floc_zo

.. currentmodule:: watertap.unit_models.zero_order.coag_and_floc_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.coag_and_floc_zo
    :members:
    :noindex:
