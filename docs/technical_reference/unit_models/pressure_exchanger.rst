Pressure Exchanger
==================
This pressure exchanger unit model:
    * is 0-dimensional
    * supports a single liquid phase only
    * supports steady-state only

Degrees of Freedom
------------------
Generally, pressure exchangers are not used in isolation and are part of an energy recovery system that includes a mixer,
splitter, and a booster pump, as shown in Figure 1. This energy recovery system only adds one degree of freedom to the overall model.

Typically, the pressure exchanger efficiency (``efficiency_pressure_exchanger``) is fixed to fully specify the energy recovery system subject to the following additional constraints:
    * volumetric flowrate is equal on both sides of the pressure exchanger (constraint is included in pressure exchanger model)
    * booster pump matches the high pressure pump outlet pressure (constraint must be added by the user or at the mixer with the equality momentum mixing type option)

Model Structure
------------------


Sets
----


Variables
----------


Equations
-----------


Class Documentation
-------------------


