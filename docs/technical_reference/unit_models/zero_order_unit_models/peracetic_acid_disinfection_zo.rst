Peracetic Acid Disinfection (ZO)
================================

Model Type
----------
This unit model is formulated as a **reactive single-inlet, double-outlet** model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

Electricity Consumption
-----------------------
Electricity consumption is calculated using the **constant_intensity** helper function.
See documentation for :ref:`Helper Methods for Electricity Demand<electricity_methods>`.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.peracetic_acid_disinfection_zo.PeraceticAcidDisinfectionData.cost_peracetic_acid` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Hydraulic retention time of water treatment unit", "HRT", ":math:`h`"
   "Average mass of an E. coli cell", "ecoli_cell_mass", ":math:`kg`"
   "Weight fraction of PAA in disinfection solution", "disinfection_solution_wt_frac_PAA", ":math:`dimensionless`"
   "Disinfection solution density", "disinfection_solution_density", ":math:`kg/l`"
   "Volumetric flowrate of disinfection solution", "disinfection_solution_flow_vol", ":math:`l/s`"
   "Volume of water treatment unit", "reactor_volume", ":math:`m^3`"
   "Concentration of E. coli at reactor inlet", "inlet_ecoli_conc", ":math:`1/l`"
   "Concentration of E. coli at reactor outlet", "outlet_ecoli_conc", ":math:`1/l`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for disinfection solution flowrate", "disinfection_solution_flow_vol_rule"
   "Constraint for reactor volume", "reactor_volume_rule"
   "Constraint for E. coli inlet concentration", "ecoli_inlet_concentration_rule"
   "Constraint for E. coli outlet concentration", "ecoli_outlet_concentration_rule"

.. index::
   pair: watertap.unit_models.zero_order.peracetic_acid_disinfection_zo;peracetic_acid_disinfection_zo

.. currentmodule:: watertap.unit_models.zero_order.peracetic_acid_disinfection_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.peracetic_acid_disinfection_zo
    :members:
    :noindex:
