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

   "Dosing rate of alum", "alum_dose", ":math:`mg/l`"
   "Dosing rate of polymer", "polymer_dose", ":math:`mg/l`"
   "Ratio of anionic to cationic polymer in dosage", "anion_to_cation_polymer_ratio", ":math:`dimensionless`"
   "Dosing rate of anionic polymer", "anionic_polymer_dose", ":math:`mg/l`"
   "Dosing rate of cationic polymer", "cationic_polymer_dose", ":math:`mg/l`"
   "Mass flow rate of chemical solution", "chemical_flow_mass", ":math:`kg/s`"
   "Mass flow rate of chemical solution", "chemical_flow_mass", ":math:`kg/s`"
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
   "Rapid Mix Power Consumption", "power_rapid_mix", ":math:`kW`"
   "Floc Power Consumption", "power_floc", ":math:`kW`"
   "Total Power Consumption", "electricity", ":math:`kW`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "None", "rapid_mix_basin_vol_constraint"
   "None", "floc_basin_vol_constraint"
   "None", "chemical_flow_constraint"
   "None", "chemical_flow_constraint"
   "None", "anionic_polymer_dose_constraint"
   "None", "cationic_polymer_dose_constraint"
   "Constraint for rapid mix power consumption", "rule_power_rapid_mix"
   "Constraint for floc power consumption", "rule_power_floc"
   "Total power consumption", "electricity_constraint"

.. index::
   pair: watertap.unit_models.zero_order.coag_and_floc_zo;coag_and_floc_zo

.. currentmodule:: watertap.unit_models.zero_order.coag_and_floc_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.coag_and_floc_zo
    :members:
    :noindex:
