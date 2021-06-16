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

Typically, the following variables are fixed, in addition to state variables at the inlet:
    * membrane water permeability, A
    * membrane salt permeability, B
    * permeate pressure
    * membrane area

On the other hand, configuring the RO unit to calculate concentration polarization effects, mass transfer
coefficient, and pressure drop would result in 7 degrees of freedom. In this case, in addition to the
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

*Solute depends on the imported property model; example shown here is for the NaCl property model.

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Solvent permeability coefficient", ":math:`A`", "A_comp", "[t, j]", ":math:`\text{m/Pa/s}`"
   "Solute permeability coefficient", ":math:`B`", "B_comp", "[t, j]", ":math:`\text{m/s}`"
   "Mass density of pure water", ":math:`\rho_w`", "dens_solvent", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Mass flux across membrane", ":math:`J`", "flux_mass_io_phase_comp", "[t, x, p, j]", ":math:`\text{kg/s}\text{/m}^2`"
   "Membrane area", ":math:`A_m`", "area", "None", ":math:`\text{m}^2`"
   "Component recovery rate", ":math:`R_j`", "recovery_mass_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Volumetric recovery rate", ":math:`R`", "recovery_vol_phase", "[t, p]", ":math:`\text{dimensionless}`"
   "Observed solute rejection", ":math:`r_s`", "rejection_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Over-pressure ratio", ":math:`P_{f,out}/Δ\pi_{out}`", "over_pressure_ratio", "[t]", ":math:`\text{dimensionless}`"

The following variables are only built when specific configuration key-value pairs are selected.

if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.fixed``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Concentration polarization modulus", ":math:`C_{f,m}/C_f`", "cp_modulus", "[t, j]", ":math:`\text{dimensionless}`"

if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Mass transfer coefficient in feed channel", ":math:`k_f`", "Kf_io", "[t, x, j]", ":math:`\text{m/s}`"

if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Feed-channel height", ":math:`h_{ch}`", "channel_height", "None", ":math:`\text{m}`"
   "Hydraulic diameter", ":math:`d_h`", "dh", "None", ":math:`\text{m}`"
   "Spacer porosity", ":math:`\epsilon_{sp}`", "spacer_porosity", "None", ":math:`\text{dimensionless}`"
   "Reynolds number", ":math:`Re`", "N_Re_io", "[t, x]", ":math:`\text{dimensionless}`"


if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Schmidt number", ":math:`Sc`", "N_Sc_io", "[t, x]", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh`", "N_Sh_io", "[t, x]", ":math:`\text{dimensionless}`"

if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is **NOT** set to ``PressureChangeType.fixed_per_stage``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane length", ":math:`L`", "length", "None", ":math:`\text{m}`"
   "Membrane width", ":math:`W`", "width", "None", ":math:`\text{m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.fixed_per_unit_length``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Pressure drop per unit length of feed channel", ":math:`ΔP/Δx`", "dP_dx", "[t]", ":math:`\text{Pa/m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Feed-channel velocity", ":math:`v_f`", "velocity_io", "[t, x]", ":math:`\text{m/s}`"
   "Friction factor", ":math:`f`", "friction_factor_darcy_io", "[t, x]", ":math:`\text{dimensionless}`"
   "Pressure drop per unit length of feed channel at inlet/outlet", ":math:`ΔP/Δx`", "dP_dx_io", "[t, x]", ":math:`\text{Pa/m}`"


Equations
-----------
#TODO:

#.. csv-table::
   #:header: "Description", "Equation"




Class Documentation
-------------------

.. automodule:: proteuslib.unit_models.reverse_osmosis_0D
   :members:
   :noindex:



