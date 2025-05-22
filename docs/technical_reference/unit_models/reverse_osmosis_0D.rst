.. _RO_0D:

Reverse Osmosis (0D)
====================

.. code-block:: python

   from watertap.unit_models.reverse_osmosis_0D import ReverseOsmosis0D

This reverse osmosis (RO) unit model
   * is 0-dimensional
   * supports a single liquid phase only
   * supports steady-state only
   * supports both solution-diffusion (SD) and Spiegler-Kedem-Katchalsky (SKK) models [1]
   * supports flat-sheet and spiral-wound module designs
   * assumes isothermal conditions

.. index::
   pair: watertap.unit_models.reverse_osmosis_0D;reverse_osmosis_0D

.. currentmodule:: watertap.unit_models.reverse_osmosis_0D

Degrees of Freedom
------------------
Aside from the inlet feed state variables (i.e. temperature, pressure, component flowrates), the RO model has
at least 4 degrees of freedom that should be fixed for the unit to be fully specified.

Typically, the following variables are fixed, in addition to state variables at the inlet:
    * membrane water permeability, A
    * membrane salt permeability, B
    * permeate pressure
    * membrane area

On the other hand, configuring the RO unit to calculate concentration polarization effects, mass transfer
coefficient, and pressure drop would result in 3 additional degrees of freedom. In this case, in addition to the
previously fixed variables, we typically fix the following variables to fully specify the unit:

    * feed-spacer porosity
    * feed-channel height
    * membrane length *or* membrane width *or* inlet Reynolds number

Model Structure
------------------
This RO model consists of 1 MembraneChannel0DBlock for the feed-side, a StateBlock indexed by time and space for the permeate-side (permeate_side[t, x]),
and a StateBlock for the final permeate at the outlet (mixed_permeate).
 
* The feed-side includes 2 StateBlocks (properties_in and properties_out) which are used for mass, energy, and momentum balances, and 2 additional StateBlocks for the conditions at the membrane interface (properties_interface_in and properties_interface_out).
* The permeate-side includes 3 StateBlocks (properties_in, properties_out, and mixed_permeate). The inlet and outlet StateBlocks are used to only determine the permeate solute concentration for solvent and solute flux at the feed-side inlet and outlet, while the mixed StateBlock is used for mass balance based on the average flux.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property model; example shown here is for the NaCl property model.

.. _0dro_variables:

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Solvent permeability coefficient", ":math:`A`", "A_comp", "[t, j]", ":math:`\text{m/Pa/s}`"
   "Solute permeability coefficient", ":math:`B`", "B_comp", "[t, j]", ":math:`\text{m/s}`"
   "Mass density of solvent", ":math:`\rho_{solvent}`", "dens_solvent", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Mass flux across membrane", ":math:`J`", "flux_mass_phase_comp", "[t, x, p, j]", ":math:`\text{kg/s}\text{/m}^2`"
   "Membrane area", ":math:`A_m`", "area", "None", ":math:`\text{m}^2`"
   "Component recovery rate", ":math:`R_j`", "recovery_mass_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Volumetric recovery rate", ":math:`R_{vol}`", "recovery_vol_phase", "[t, p]", ":math:`\text{dimensionless}`"
   "Observed solute rejection", ":math:`r_j`", "rejection_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Over-pressure ratio", ":math:`P_{f,out}/Δ\pi_{out}`", "over_pressure_ratio", "[t]", ":math:`\text{dimensionless}`"
   "Mass transfer to permeate", ":math:`M_p`", "mass_transfer_phase_comp", "[t, p, j]", ":math:`\text{kg/s}`"
   "Feed-channel volume", ":math:`V_{f}`", "feed_side.volume", "[t]", ":math:`\text{m}^3`"

The following variables are only built when specific configuration key-value pairs are selected.

if ``has_pressure_change`` is set to ``True``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Pressure drop", ":math:`ΔP`", "deltaP", "[t]", ":math:`\text{Pa}`"

if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.fixed``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Concentration polarization modulus", ":math:`CP_{mod}`", "feed_side.cp_modulus", "[t, j]", ":math:`\text{dimensionless}`"

if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Mass transfer coefficient in feed channel", ":math:`k_f`", "feed_side.K", "[t, x, j]", ":math:`\text{m/s}`"

if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Feed-channel height", ":math:`h_{ch}`", "feed_side.channel_height", "None", ":math:`\text{m}`"
   "Hydraulic diameter", ":math:`d_h`", "feed_side.dh", "None", ":math:`\text{m}`"
   "Spacer porosity", ":math:`\epsilon_{sp}`", "feed_side.spacer_porosity", "None", ":math:`\text{dimensionless}`"
   "Reynolds number", ":math:`Re`", "feed_side.N_Re", "[t, x]", ":math:`\text{dimensionless}`"


if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Schmidt number", ":math:`Sc`", "feed_side.N_Sc_comp", "[t, x, j]", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh`", "feed_side.N_Sh_comp", "[t, x, j]", ":math:`\text{dimensionless}`"

if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is **NOT** set to ``PressureChangeType.fixed_per_stage``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane length", ":math:`L`", "length", "None", ":math:`\text{m}`"
   "Membrane width", ":math:`W`", "width", "None", ":math:`\text{m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.fixed_per_unit_length``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Average pressure drop per unit length of feed channel", ":math:`(\frac{ΔP}{Δx})_{avg}`", "feed_side.dP_dx", "[t]", ":math:`\text{Pa/m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Feed-channel velocity", ":math:`v_f`", "feed_side.velocity", "[t, x]", ":math:`\text{m/s}`"
   "Friction factor", ":math:`f`", "feed_side.friction_factor_darcy", "[t, x]", ":math:`\text{dimensionless}`"
   "Pressure drop per unit length of feed channel at inlet/outlet", ":math:`ΔP/Δx`", "feed_side.dP_dx", "[t, x]", ":math:`\text{Pa/m}`"

if ``transport_model`` is set to ``TransportModel.SKK``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Reflection coefficient", ":math:`\sigma`", "reflect_coeff", "None", ":math:`\text{dimensionless}`"
   "Alpha", ":math:`\alpha`", "alpha", "None", ":math:`\text{s/m}`"

.. _0dro_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "Solvent flux across membrane (solution-diffusion)", ":math:`J_{solvent} = \rho_{solvent} A(P_{f} - P_p - (\pi_{f}-\pi_{p}))`"
   "Solvent flux across membrane (SKK)", ":math:`J_{solvent} = \rho_{solvent} A(P_{f} - P_p - \sigma(\pi_{f}-\pi_{p}))`"
   "Solute flux across membrane (solution-diffusion)", ":math:`J_{solute} = B(C_{f} - C_{p})`"
   "Solute flux across membrane (SKK)", ":math:`J_{solute} = B(C_{f} - C_{p}) + (1 - \sigma)\frac{J_{solvent}}{\rho_{solvent}}C_{f}`"
   "Alpha (SKK)", ":math:`\alpha = \frac{1 - \sigma}{B}`" 
   "Average flux across membrane", ":math:`J_{avg, j} = \frac{1}{2}\sum_{x} J_{x, j}`"
   "Permeate mass flow by component j", ":math:`M_{p, j} = A_m J_{avg,j}`"
   "Permeate-side solute mass fraction", ":math:`X_{x, j} = \frac{J_{x, j}}{\sum_{x} J_{x, j}}`"
   "Feed-side membrane-interface solute concentration", ":math:`C_{interface} = CP_{mod}C_{bulk}=C_{bulk}\exp(\frac{J_{solvent}}{k_f})-\frac{J_{solute}}{J_{solvent}}(\exp(\frac{J_{solvent}}{k_f})-1)`"
   "Concentration polarization modulus",":math:`CP_{mod} = C_{interface}/C_{bulk}`"
   "Mass transfer coefficient",":math:`k_f = \frac{D Sh}{d_h}`"
   "Sherwood number",":math:`Sh = 0.46 (Re Sc)^{0.36}`"
   "Schmidt number",":math:`Sc = \frac{\mu}{\rho D}`"
   "Reynolds number",":math:`Re = \frac{\rho v_f d_h}{\mu}`"
   "Hydraulic diameter",":math:`d_h = \frac{4\epsilon_{sp}}{2/h_{ch} + (1-\epsilon_{sp})8/h_{ch}}`"
   "Cross-sectional area",":math:`A_c = h_{ch}W\epsilon_{sp}`"
   "Membrane area (flat-plate)",":math:`A_m = LW`"
   "Membrane area (spiral-wound)",":math:`A_m = 2LW`"
   "Pressure drop",":math:`ΔP = (\frac{ΔP}{Δx})_{avg}L`"
   "Feed-channel velocity",":math:`v_f = Q_f/A_c`"
   "Feed-channel volume", ":math:`V_{f} = LWh_{ch}`"
   "Friction factor (flat-plate)",":math:`f = 0.42+\frac{189.3}{Re}`"
   "Friction factor (spiral-wound)",":math:`f = \frac{6.23}{Re^{0.3}}`"
   "Pressure drop per unit length",":math:`\frac{ΔP}{Δx} = \frac{1}{2d_h}f\rho v_f^{2}`"
   "Component recovery rate",":math:`R_j = \frac{M_{p,j}}{M_{f,in,j}}`"
   "Volumetric recovery rate",":math:`R_{vol} = \frac{Q_{p}}{Q_{f,in}}`"
   "Observed solute rejection", ":math:`r_j = 1 - \frac{C_{p,mix}}{C_{f,in}}`" 

Class Documentation
-------------------

* :mod:`watertap.unit_models.reverse_osmosis_0D`

References
----------
Spiegler, K. S., & Kedem, O. (1966). Thermodynamics of hyperfiltration (reverse osmosis): criteria for efficient membranes. Desalination, 1(4), 311-326.