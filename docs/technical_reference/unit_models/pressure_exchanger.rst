Pressure Exchanger
==================
This pressure exchanger unit model:
    * is 0-dimensional
    * is isothermal
    * supports a single liquid phase only
    * supports steady-state only
    * assumes no mixing or leakage between the low and high pressure side
    * assumes equal flowrates on both sides

.. index::
   pair: watertap.unit_models.pressure_exchanger;pressure_exchanger

.. currentmodule:: watertap.unit_models.pressure_exchanger

Degrees of Freedom
------------------
Generally, pressure exchangers are not used in isolation and are part of an energy recovery system that includes a mixer,
splitter, and a booster pump, as shown in Figure 1. This energy recovery system only adds one degree of freedom to the overall model.

Typically, the following variables are fixed to fully specify the energy recovery system:
    * the pressure exchanger efficiency (`efficiency_pressure_exchanger` located on the pressure exchanger unit model)
    * the booster pump efficiency (`efficiency_pump` located on the pump unit model)

Where the system is also subject to following constraints:
    * volumetric flowrate is equal on both sides of the pressure exchanger (constraint is included in pressure exchanger unit model)
    * booster pump matches the high pressure pump outlet pressure (constraint must be added by the user or at the mixer with the equality momentum mixing type option)

.. figure:: ../../_static/unit_models/energy_recovery_system.png
    :width: 600
    :align: center
    
    Figure 1. Schematic representation of an energy recovery system using a pressure exchanger.

Model Structure
------------------
The pressure exchanger model consists of 2 `ControlVolume0DBlocks`: one for the low-pressure side and high-pressure side.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property package; example shown here is for the NaCl property model.

Variables
----------

The pressure exchanger unit model includes the following variables:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Pyomo Type" 

   "Efficiency", ":math:`\eta`", "efficiency_pressure_exchanger", "[t]", ":math:`\text{dimensionless}`", "Var"

Each control volume (i.e. `low_presssure_side`, and `high_pressure_side`) has the following variables of interest:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Pyomo Type" 

   "Pressure change", ":math:`ΔP`", "deltaP", "[t]", "\*", "Var"
   "Work", ":math:`W`", "work", "[t]", "\*", "Expression"

\*Units depends on the imported property package

Each property block on both control volumes (i.e. `properties_in` and `properties_out`) has the following variables of interest:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Pyomo Type" 

   "Volumetric flowrate", ":math:`Q`", "flow_vol", "None", "\*", "Var"
   "Temperature", ":math:`T`", "temperature", "[t]", "\*", "Var"
   "Pressure", ":math:`P`", "pressure", "[t]", "\*", "Var"

\*Units depends on the imported property package
   

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "Mass balance for each side", ":math:`M_{out, j} = M_{in, j}`"
   "Momentum balance for each side", ":math:`P_{out} = P_{in} + ΔP`"
   "Isothermal assumption for each side", ":math:`T_{out} = T_{in}`"
   "Equal volumetric flowrate*", ":math:`Q_{LPS} = Q_{HPS}`"
   "Pressure transfer*", ":math:`ΔP_{LPS} = - \eta ΔP_{HPS}`"

\* LPS stands for low pressure side, HPS stands for high pressure side

Class Documentation
-------------------

* :mod:`watertap.unit_models.pressure_exchanger`
