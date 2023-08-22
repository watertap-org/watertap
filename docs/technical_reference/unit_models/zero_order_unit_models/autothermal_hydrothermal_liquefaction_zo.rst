Autothermal Hydrothermal Liquefaction (ZO)
==========================================

Model Type
----------
This unit model is formulated as a **reactive single-inlet, double-outlet** model form.
See documentation for :ref:`reactive single-inlet, double-outlet Helper Methods<sidor_methods>`.

Electricity Consumption
-----------------------
The constraint used to calculate energy consumption is described in the Additional Constraints section below. More details can be found in the unit model class.

Costing Method
--------------
Costing is calculated using the :py:meth:`~watertap.unit_models.zero_order.autothermal_hydrothermal_liquefaction_zo.ATHTLZOData.cost_autothermal_hydrothermal_liquefaction` method.
For full details on costing, see documentation for the :ref:`zero-order costing package<zero_order_costing>`.

Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Variable Name", "Units"

   "Inlet mass flowrate", "flow_mass_in", ":math:`t/h`"
   "Electricity consumption of unit", "electricity", ":math:`kW`"
   "Electricity intensity with respect to inlet flowrate", "energy_electric_flow_mass", ":math:`kWh/t`"
   "Dosage of catalyst per inlet flow", "catalyst_dosage", ":math:`lb/t`"
   "Catalyst flow", "catalyst_flow", ":math:`lb/h`"

Additional Constraints
----------------------

.. csv-table::
   :header: "Description", "Constraint Name"

   "Constraint for inlet mass flowrate.", "cons_flow_mass"
   "Constraint for electricity consumption based on inlet flow rate.", "electricity_consumption"
   "Constraint for catalyst flow based on inlet flow rate.", "eq_catalyst_flow"

.. index::
   pair: watertap.unit_models.zero_order.autothermal_hydrothermal_liquefaction_zo;autothermal_hydrothermal_liquefaction_zo

.. currentmodule:: watertap.unit_models.zero_order.autothermal_hydrothermal_liquefaction_zo

Class Documentation
-------------------

.. automodule:: watertap.unit_models.zero_order.autothermal_hydrothermal_liquefaction_zo
    :members:
    :noindex:
