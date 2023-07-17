Modular Encapsulated Two-stage Anaerobic Biological Reactor (ZO)
================================================================

Model Type
----------
This unit model is formulated as a **reactive single-inlet, double-outlet** model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.metab_zo.MetabZOData.cost_metab` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Reactor volume", "volume", ":math:`m^3`"
   "Hydraulic residence time", "hydraulic_retention_time", ":math:`h`"
   "Electricity demand of unit", "electricity", ":math:`kW`"
   "Thermal demand of unit", "heat", ":math:`kW`"
   "Electricity intensity of mixer with respect to reactor volume", "energy_electric_mixer_vol", ":math:`kW/m^3`"
   "Electricity intensity of vacuum pump with respect to product gas flow", "energy_electric_vacuum_flow_vol_byproduct", ":math:`h*kW/kg`"
   "Thermal energy intensity of reactor with respect to inlet volumetric flowrate", "energy_thermal_flow_vol_inlet", ":math:`kJ/m^3`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for reactor volume based on hydraulic residence time", "eq_reactor_volume"
   "Constraint for electricity consumption based on feed flowrate.", "electricity_consumption"
   "Constraint for heat demand based on feed flowrate.", "heat_demand"

.. index::
   pair: watertap.unit_models.zero_order.metab_zo;metab_zo

.. currentmodule:: watertap.unit_models.zero_order.metab_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.metab_zo
    :members:
    :noindex:
