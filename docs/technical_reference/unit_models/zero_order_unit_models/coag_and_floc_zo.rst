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
Costing is calculated using the **cost_coag_and_floc** method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Dosing rate of alum", "alum_dose", "mg/l"
   "Dosing rate of polymer", "polymer_dose", "mg/l"
   "Ratio of anionic to cationic polymer in dosage", "anion_to_cation_polymer_ratio", "dimensionless"
   "Dosing rate of anionic polymer", "anionic_polymer_dose", "mg/l"
   "Dosing rate of cationic polymer", "cationic_polymer_dose", "mg/l"
   "Mass flow rate of chemical solution", "chemical_flow_mass", "kg/s"
   "Mass flow rate of chemical solution", "chemical_flow_mass", "kg/s"
   "Rapid Mix Retention Time", "rapid_mix_retention_time", "s"
   "Floc Retention Time", "floc_retention_time", "min"
   "Rapid Mix Basin Volume", "rapid_mix_basin_vol", "m**3"
   "Floc Basin Volume", "floc_basin_vol", "m**3"
   "Number of Rapid Mixers", "num_rapid_mixers", "dimensionless"
   "Number of Floc Mixers", "num_floc_mixers", "dimensionless"
   "Number of Rapid Mix Processes", "num_rapid_mix_processes", "dimensionless"
   "Number of Floc Processes", "num_floc_processes", "dimensionless"
   "Number of Coagulation Processes", "num_coag_processes", "dimensionless"
   "Number of Floc Injection Processes", "num_floc_injection_processes", "dimensionless"
   "Rapid Mix Velocity Gradient", "velocity_gradient_rapid_mix", "1/s"
   "Floc Velocity Gradient", "velocity_gradient_floc", "1/s"
   "Rapid Mix Power Consumption", "power_rapid_mix", "kW"
   "Floc Power Consumption", "power_floc", "kW"
   "Total Power Consumption", "electricity", "kW"

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
