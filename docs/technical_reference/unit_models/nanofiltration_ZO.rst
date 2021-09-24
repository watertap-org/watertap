Nanofiltration (ZO)
====================
This nanofiltration (NF) unit model
   * is a zero-order model that enables the user to specify performance in terms of membrane solvent flux and solute/ion rejection
   * supports a single liquid phase only
   * supports steady-state only
   * assumes isothermal conditions

.. index::
   pair: proteuslib.unit_models.nanofiltration_ZO;nanofiltration_ZO

.. currentmodule:: proteuslib.unit_models.nanofiltration_ZO

Degrees of Freedom
------------------
The zero-order NF model has at least 7 degrees of freedom that should be fixed for the unit to be fully specified.
Typically, the following variables are fixed:

   * solvent mass flow rate at the inlet
   * inlet temperature
   * inlet pressure
   * solvent volumetric flux of the membrane
   * permeate pressure

There are 2 degrees of freedom for each solute in a given property model:

   * solute rejection of the membrane
   * solute mass flow rate at the inlet

Note, when a set of solutes comprises ions and includes more than one ion pair, the solute rejection of one ion could be left unfixed
to satisfy an electroneutrality constraint. In this case, the last degree of freedom can be eliminated by fixing membrane area.

Model Structure
------------------
This NF model consists of 1 ControlVolume0DBlock for the feed-side of the membrane and 1 StateBlock (properties_permeate) for the permeate exiting the NF membrane module.
The feed-side includes 2 StateBlocks (properties_in and properties_out) which are used for mass and momentum balances.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'Solute']*"

\*Solute depends on the imported property model, and more than one solute can be provided; ions can be represented as solutes in a given property model to be recognized by this NF model.

.. _ZO_NF_variables:

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Mass density of pure water", ":math:`\rho_w`", "dens_solvent", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Solvent volumetric flux across membrane", ":math:`J_{solv}`", "flux_vol_solvent", "[t, j]", ":math:`\text{m}^3\text{/m}^2\text{/s}`"
   "Membrane area", ":math:`A_m`", "area", "None", ":math:`\text{m}^2`"
   "Component recovery rate", ":math:`R_j`", "recovery_mass_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Volumetric recovery rate", ":math:`R_{vol}`", "recovery_vol_phase", "[t, p]", ":math:`\text{dimensionless}`"
   "Observed solute rejection", ":math:`r_j`", "rejection_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Mass transfer to permeate", ":math:`M_p`", "mass_transfer_phase_comp", "[t, p, j]", ":math:`\text{kg/s}`"

.. _ZO_NF_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "Solvent flux across membrane", ":math:`J_{solvent} = \rho_{solvent} A(P_{f} - P_p - (\pi_{f}-\pi_{p}))`"
   "Solute flux across membrane", ":math:`J_{solute} = B(C_{f} - C_{p})`"
   "Average flux across membrane", ":math:`J_{avg, j} = \frac{1}{2}\sum_{x} J_{x, j}`"
   "Permeate mass flow by component j", ":math:`M_{p, j} = A_m J_{avg,j}`"
   "Permeate-side membrane-interface solute mass fraction", ":math:`X_{x, j} = \frac{J_{x, j}}{\sum_{x} J_{x, j}}`"
   "Feed-side membrane-interface solute concentration", ":math:`C_{interface} = CP_{mod}C_{bulk}=C_{bulk}\exp(\frac{J_{solvent}}{k_f})-\frac{J_{solute}}{J_{solvent}}(\exp(\frac{J_{solvent}}{k_f})-1)`"
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

.. automodule:: proteuslib.unit_models.nanofiltration_ZO
    :members:
    :noindex:

