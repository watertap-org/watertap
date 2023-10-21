Electrocoagulation (ZO)
=======================

Model Type
----------
This unit model is formulated as a **single-input, double-output** model form.
See documentation for :ref:`single-input, double-output Helper Methods<sido_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.electrocoagulation_zo.ElectrocoagulationZOData.cost_electrocoagulation` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Cathode area", "cathode_area", ":math:`m^2`"
   "Anode area", "anode_area", ":math:`m^2`"
   "Electrode thickness", "electrode_thick", ":math:`m`"
   "Electrode mass", "electrode_mass", ":math:`kg`"
   "Electrode volume", "electrode_volume", ":math:`m^3`"
   "Electrode gap", "electrode_gap", ":math:`m`"
   "Conductivity", "conductivity", ":math:`S/m`"
   "Applied current", "applied_current", ":math:`A`"
   "Current efficiency", "current_efficiency", ":math:`dimensionless`"
   "Cell voltage", "cell_voltage", ":math:`V`"
   "Overpotential", "overpotential", ":math:`V`"
   "Reactor volume", "reactor_volume", ":math:`m^3`"
   "Metal dose", "metal_dose", ":math:`kg/L`"
   "Ohmic resistance of solution", "ohmic_resistance", ":math:`Ω`"
   "Charge loading rate", "charge_loading_rate", ":math:`C/l`"
   "Current density", "current_density", ":math:`A/m^2`"
   "Power required", "power_required", ":math:`W`"
   "Floc basin volume", "floc_basin_vol", ":math:`m^3`"
   "Floc retention time", "floc_retention_time", ":math:`min`"
   "Overpotential calculation", "eq_overpotential"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Charge loading rate equation", "eq_charge_loading_rate"
   "Total current required", "eq_applied_current"
   "Total electrode area required", "eq_electrode_area_total"
   "Cell voltage", "eq_cell_voltage"
   "Electrode volume", "eq_electrode_volume"
   "Cathode/Anode area", "eq_cathode_anode"
   "Total reactor volume", "eq_reactor_volume"
   "Total flocculation tank volume", "eq_floc_reactor_volume"
   "Ohmic resistance", "eq_ohmic_resistance"
   "Electrode mass", "eq_electrode_mass"
   "Power required", "eq_power_required"

.. index::
   pair: watertap.unit_models.zero_order.electrocoagulation_zo;electrocoagulation_zo

.. currentmodule:: watertap.unit_models.zero_order.electrocoagulation_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.electrocoagulation_zo
    :members:
    :noindex:
