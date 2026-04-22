.. _OARO_0D:

Osmotically Assisted Reverse Osmosis (0D)
=========================================

.. code-block:: python

   from watertap.unit_models.osomtically_assisted_reverse_osmosis_0D import OsmoticallyAssistedReverseOsmosis0D

This osmotically assisted reverse osmosis (OARO) unit model
   * is 0-dimensional
   * supports a single liquid phase only
   * supports steady-state only
   * supports flat-sheet module designs
   * is based on the solution-diffusion model and film theory
   * assumes isothermal conditions
   * assumes the feed-side flows in the forward direction
   * assumes the permeate-side flows in the backwards direction

.. index::
   pair: watertap.unit_models.osmotically_assisted_reverse_osmosis_0D;osomtically_assisted_reverse_osmosis_0D

.. currentmodule:: watertap.unit_models.osmotically_assisted_reverse_osmosis_0D

Degrees of Freedom
------------------
Aside from the feed-side and permeate-side inlet state variables (i.e. temperature, pressure, component flowrates), the OARO model has
at least 4 degrees of freedom that should be fixed for the unit to be fully specified. Unlike RO, which only
accounts for concentration polarization on the feed side, the OARO model includes a structural parameter
variable, which is used to calculate the membrane-interface concentration on the permeate side.

Typically, the following variables are fixed for the OARO model, in addition to state variables at the inlet:
    * membrane water permeability, A
    * membrane salt permeability, B
    * membrane area
    * structural parameter

On the other hand, configuring the OARO unit to calculate concentration polarization effects, mass transfer
coefficient, and pressure drop would result in 6 additional degrees of freedom. In this case, in addition to the
previously fixed variables, we typically fix the following variables to fully specify the unit:

    * feed-side spacer porosity
    * feed-side channel height
    * permeate-side space porosity
    * permeate-side channel height
    * feed-side Reynolds number *or* water mass recovery

Model Structure
------------------
This OARO model consists of a separate MembraneChannel0DBlock for the feed-side and the permeate-side of the membrane.

* The feed-side includes 2 StateBlocks (properties_in and properties_out) which are used for mass, energy, and momentum balances, and 2 additional StateBlocks for the conditions at the membrane interface (properties_interface_in and properties_interface_out).
* The permeate-side includes 2 StateBlocks (properties_in and properties_out) which are used for mass, energy, and momentum balances, and 2 additional StateBlocks for the conditions at the membrane interface (properties_interface_in and properties_interface_out).

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property model; example shown here is for the NaCl property model.

.. _0doaro_variables:

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
   "Membrane structural parameter", ":math:`S`", "structural_parameter", "[None]", ":math:`\text{m}`"
   "Mass transfer to permeate", ":math:`M_p`", "mass_transfer_phase_comp", "[t, p, j]", ":math:`\text{kg/s}`"


The following variables are only built when specific configuration key-value pairs are selected.

if ``has_pressure_change`` is set to ``True``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Pressure drop", ":math:`ΔP`", "deltaP", "[t]", ":math:`\text{Pa}`"

if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.fixed``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Concentration polarization modulus in feed channel", ":math:`CP_{mod,f}`", "feed_side.cp_modulus", "[t, j]", ":math:`\text{dimensionless}`"
   "Concentration polarization modulus in permeate channel", ":math:`CP_{mod,p}`", "permeate_side.cp_modulus", "[t, j]", ":math:`\text{dimensionless}`"

if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Mass transfer coefficient in feed channel", ":math:`k_f`", "feed_side.K", "[t, x, j]", ":math:`\text{m/s}`"
   "Mass transfer coefficient in permeate channel", ":math:`k_p`", "permeate_side.K", "[t, x, j]", ":math:`\text{m/s}`"

if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Feed-channel height", ":math:`h_{ch,f}`", "feed_side.channel_height", "None", ":math:`\text{m}`"
   "Hydraulic diameter", ":math:`d_{h,f}`", "feed_side.dh", "None", ":math:`\text{m}`"
   "Spacer porosity", ":math:`\epsilon_{sp,f}`", "feed_side.spacer_porosity", "None", ":math:`\text{dimensionless}`"
   "Reynolds number", ":math:`Re_{f}`", "feed_side.N_Re", "[t, x]", ":math:`\text{dimensionless}`"
   "Permeate-channel height", ":math:`h_{ch,p}`", "permeate_side.channel_height", "None", ":math:`\text{m}`"
   "Hydraulic diameter", ":math:`d_{h,p}`", "permeate_side.dh", "None", ":math:`\text{m}`"
   "Spacer porosity", ":math:`\epsilon_{sp,p}`", "permeate_side.spacer_porosity", "None", ":math:`\text{dimensionless}`"
   "Reynolds number", ":math:`Re_{p}`", "permeate_side.N_Re", "[t, x]", ":math:`\text{dimensionless}`"


if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Schmidt number", ":math:`Sc_f`", "feed_side.N_Sc", "[t, x]", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh_f`", "feed_side.N_Sh", "[t, x]", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc_p`", "permeate_side.N_Sc", "[t, x]", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh_p`", "permeate_side.N_Sh", "[t, x]", ":math:`\text{dimensionless}`"

if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is **NOT** set to ``PressureChangeType.fixed_per_stage``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane length", ":math:`L`", "length", "None", ":math:`\text{m}`"
   "Membrane width", ":math:`W`", "width", "None", ":math:`\text{m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.fixed_per_unit_length``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Average pressure drop per unit length of feed channel", ":math:`(\frac{ΔP}{Δx})_{avg,f}`", "feed_side.dP_dx", "[t]", ":math:`\text{Pa/m}`"
   "Average pressure drop per unit length of permate channel", ":math:`(\frac{ΔP}{Δx})_{avg,p}`", "permeate_side.dP_dx", "[t]", ":math:`\text{Pa/m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Feed-channel velocity", ":math:`v_f`", "feed_side.velocity", "[t, x]", ":math:`\text{m/s}`"
   "Friction factor", ":math:`f`", "feed_side.friction_factor_darcy", "[t, x]", ":math:`\text{dimensionless}`"
   "Pressure drop per unit length of feed channel at inlet/outlet", ":math:`(ΔP/Δx)_f`", "feed_side.dP_dx", "[t, x]", ":math:`\text{Pa/m}`"
   "Permeate-channel velocity", ":math:`v_p`", "permeate_side.velocity", "[t, x]", ":math:`\text{m/s}`"
   "Pressure drop per unit length of permeate channel at inlet/outlet", ":math:`(ΔP/Δx)_p`", "permeate_side.dP_dx", "[t, x]", ":math:`\text{Pa/m}`"

.. _0doaro_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "Solvent flux across membrane", ":math:`J_{solvent} = \rho_{solvent} A(P_{f} - P_p - (\pi_{f}-\pi_{p}))`"
   "Solute flux across membrane", ":math:`J_{solute} = B(C_{f} - C_{p})`"
   "Average flux across membrane", ":math:`J_{avg, j} = \frac{1}{2}\sum_{x} J_{x, j}`"
   "Permeate mass flow by component j", ":math:`M_{p, j} = A_m J_{avg,j}`"
   "Permeate-side solute mass fraction", ":math:`X_{x, j} = \frac{J_{x, j}}{\sum_{x} J_{x, j}}`"
   "Feed-side membrane-interface solute concentration", ":math:`C_{interface,f} = CP_{mod,f}C_{bulk,f}=C_{bulk,f}\exp(\frac{J_{solvent}}{k_f})-\frac{J_{solute}}{J_{solvent}}(\exp(\frac{J_{solvent}}{k_f})-1)`"
   "Permeate-side membrane-interface solute concentration", ":math:`C_{interface,p} = CP_{mod,p}C_{bulk,p}=C_{bulk,p}\exp(-J_{solvent}(\frac{S}{D} + \frac{1}{k_p})) + \frac{J_{solute}}{J{solvent}}(1-\exp(-J_{solvent}(\frac{S}{D} + \frac{1}{k_p})))`"
   "Concentration polarization modulus",":math:`CP_{mod} = C_{interface}/C_{bulk}`"
   "Mass transfer coefficient",":math:`k_f = \frac{D Sh}{d_h}`"
   "Sherwood number",":math:`Sh = 0.46 (Re Sc)^{0.36}`"
   "Schmidt number",":math:`Sc = \frac{\mu}{\rho D}`"
   "Reynolds number",":math:`Re = \frac{\rho v_f d_h}{\mu}`"
   "Hydraulic diameter",":math:`d_h = \frac{4\epsilon_{sp}}{2/h_{ch} + (1-\epsilon_{sp})8/h_{ch}}`"
   "Cross-sectional area",":math:`A_c = h_{ch}W\epsilon_{sp}`"
   "Membrane area",":math:`A_m = LW`"
   "Pressure drop",":math:`ΔP = (\frac{ΔP}{Δx})_{avg}L`"
   "Feed-channel velocity",":math:`v_f = Q_f/A_c`"
   "Friction factor",":math:`f = 0.42+\frac{189.3}{Re}`"
   "Pressure drop per unit length",":math:`\frac{ΔP}{Δx} = \frac{1}{2d_h}f\rho v_f^{2}`"
   "Component recovery rate",":math:`R_j = \frac{M_{p,j}}{M_{f,in,j}}`"
   "Volumetric recovery rate",":math:`R_{vol} = \frac{Q_{p}}{Q_{f,in}}`"
   "Observed solute rejection", ":math:`r_j = 1 - \frac{C_{p,mix}}{C_{f,in}}`"



Class Documentation
-------------------

* :mod:`watertap.unit_models.osmotically_assisted_reverse_osmosis_0D`
