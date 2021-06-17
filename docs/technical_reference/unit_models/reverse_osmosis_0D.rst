Reverse Osmosis Unit (0D)
=========================
This reverse osmosis (RO) unit model
   * is 0-dimensional
   * supports a single liquid phase only
   * supports steady-state only
   * is based on the solution-diffusion model and film theory

.. index::
   pair: proteuslib.unit_models.reverse_osmosis_0D;reverse_osmosis_0D

.. currentmodule:: proteuslib.unit_models.reverse_osmosis_0D

Degrees of Freedom
------------------
Aside from the feed temperature, feed pressure, and component mass flow rates at the inlet, the RO model typically has
at least 4 degrees of freedom that should be fixed for the unit to be fully specified.

In Example A, the following variables were fixed, in addition to state variables at the inlet:
    * membrane water permeability, A
    * membrane salt permeability, B
    * permeate pressure
    * membrane area

The degrees of freedom will depend on which RO configuration options are selected. For example, setting
``has_pressure_change= True`` adds 1 degree of freedom. In this case, the pressure drop ``deltaP`` would need to be fixed
to eliminate that degree of freedom.

On the other hand, in Example B, configuring the RO unit to calculate concentration polarization effects, mass transfer
coefficient, and pressure drop would result in 3 more degrees of freedom than Example A. In this case, in addition to the
previously fixed variables, we typically fix the following variables to fully specify the unit:

    * feed-spacer porosity
    * feed-channel height
    * membrane length *or* membrane width

Model Structure
------------------
This RO model consists of a ControlVolume0DBlock for the feed-channel of the membrane and a StateBlock for the permeate channel.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

*Solute depends on the imported property model; example shown here is for the NaCl property model.*


Variables
----------


.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Solvent permeability coefficient", ":math:`A`", "A_comp", "[t, j]", ":math:`\text{m/Pa/s}`"
   "Solvent permeability coefficient", ":math:`B`", "B_comp", "[t, j]", ":math:`\text{m/s}`"
   "Mass density of pure water", ":math:`\rho_w`", "dens_solvent", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Mass flux across membrane", ":math:`J`", "flux_mass_io_phase_comp", "[t, x, p, j]", ":math:`\text{kg/s}\text{/m}^2`"
   "Membrane area", ":math:`A_m`", "area", "None", ":math:`\text{m}^2`"
   "Component recovery rate", ":math:`R_j`", "recovery_mass_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Volumetric recovery rate", ":math:`R`", "recovery_vol_phase", "[t, p]", ":math:`\text{dimensionless}`"


Constraints
-----------
#TODO:

.. .. csv-table::
..   :header: "Description", "Equation"




Class Documentation
-------------------


.. autoclass:: ReverseOsmosis0D
   :members:


.. autoclass:: ReverseOsmosisData
   :members:




