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

   "Electrode width", "electrode_width", ":math:`m`"
   "Electrode height", "electrode_height", ":math:`m`"
   "Electrode thickness", "electrode_thick", ":math:`m`"
   "Electrode mass", "electrode_mass", ":math:`kg`"
   "Total electrode area", "electrode_area_total", ":math:`m^2`"
   "Electrode area", "electrode_area_per", ":math:`m^2`"
   "Electrode volume", "electrode_volume_per", ":math:`m^3`"
   "Electrode gap", "electrode_gap", ":math:`m`"
   "Electrolysis time", "electrolysis_time", ":math:`min`"
   "Number of electrode pairs", "number_electrode_pairs", ":math:`dimensionless`"
   "Number of cells", "number_cells", ":math:`dimensionless`"
   "Applied current", "applied_current", ":math:`A`"
   "Current efficiency", "current_efficiency", ":math:`dimensionless`"
   "Cell voltage", "cell_voltage", ":math:`V`"
   "Overpotential", "overpotential", ":math:`V`"
   "Reactor volume total (electrochemical + flotation + sedimentation)", "reactor_volume", ":math:`m^3`"
   "Metal loading", "metal_loading", ":math:`kg/l`"
   "Ohmic resistance of solution", "ohmic_resistance", ":math:`Ω`"
   "Charge loading rate", "charge_loading_rate", ":math:`C/l`"
   "Current density", "current_density", ":math:`A/m^2`"
   "Power required", "power_required", ":math:`W`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Charge loading rate equation", "eq_charge_loading_rate"
   "Metal loading equation", "eq_metal_loading_rate"
   "Total current required", "eq_applied_current"
   "Total electrode area required", "eq_electrode_area_total"
   "Cell voltage", "eq_cell_voltage"
   "Area per electrode", "eq_electrode_area_per"
   "Electrode width", "eq_electrode_width"
   "Electrode height", "eq_electrode_height"
   "Electrode volume", "eq_electrode_volume_per"
   "Total reactor volume", "eq_reactor_volume"
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
