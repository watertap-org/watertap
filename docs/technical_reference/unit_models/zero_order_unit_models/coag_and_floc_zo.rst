Coagulation Flocculation (ZO)
=============================

Model Type
----------
This unit model is formulated as a pass-through model form.
See documentation for :ref:`pass-through Helper Methods<pt_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the f(x) helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the cost_coag_and_floc method in the zero-order costing package.
See documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

..csv-table::
   :header: "Description", "Variable Name"
   "Dosing rate of alum", "alum_dose"
   "Dosing rate of polymer", "polymer_dose"
   "Ratio of anionic to cationic polymer in dosage", "anion_to_cation_polymer_ratio"
   "Dosing rate of anionic polymer", "anionic_polymer_dose"
   "Dosing rate of cationic polymer", "cationic_polymer_dose"
   "Mass flow rate of chemical solution", "chemical_flow_mass"
   "Mass flow rate of chemical solution", "chemical_flow_mass"
   "Rapid Mix Retention Time", "rapid_mix_retention_time"
   "Floc Retention Time", "floc_retention_time"
   "Rapid Mix Basin Volume", "rapid_mix_basin_vol"
   "Floc Basin Volume", "floc_basin_vol"
   "Number of Rapid Mixers", "num_rapid_mixers"
   "Number of Floc Mixers", "num_floc_mixers"
   "Number of Rapid Mix Processes", "num_rapid_mix_processes"
   "Number of Floc Processes", "num_floc_processes"
   "Number of Coagulation Processes", "num_coag_processes"
   "Number of Floc Injection Processes", "num_floc_injection_processes"
   "Rapid Mix Velocity Gradient", "velocity_gradient_rapid_mix"
   "Floc Velocity Gradient", "velocity_gradient_floc"
   "Rapid Mix Power Consumption", "power_rapid_mix"
   "Floc Power Consumption", "power_floc"
   "Total Power Consumption", "electricity"

.. index::
   pair: watertap.unit_models.zero_order.coag_and_floc_zo;coag_and_floc_zo

.. currentmodule:: watertap.unit_models.zero_order.coag_and_floc_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.coag_and_floc_zo
    :members:
    :noindex:
