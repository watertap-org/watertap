.. _nanofiltration_DSPMDE:

Nanofiltration- Donnan Steric Pore Model with Dielectric Exclusion (0D)
=======================================================================

.. code-block:: python

   from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D

This unit model implements the Donnan Steric Pore Model with Dielectric Exclusion (DSPM-DE) for nanofiltration.

.. note::

    Documentation for the DSPM-DE model is undergoing refinement.

Model Structure
------------------
The model consists of 1 ControlVolume0DBlock for the feed-side of the membrane and includes 11 StateBlocks overall.

* The feed-side includes 2 StateBlocks (properties_in and properties_out) which are used for mass, energy, and momentum balances, and 2 additional StateBlocks for the conditions at the membrane interface (properties_interface).
* 2 StateBlocks are attributed to the membrane pore entrance (pore_entrance) at the inlet and outlet of the module.
* 2 StateBlocks are attributed to the membrane pore exit (pore_exit) at the inlet and outlet of the module.
* The permeate-side of the membrane includes 3 StateBlocks; 2 of them are attributed to the permeate side at the inlet and outlet (permeate_side), and the 3rd StateBlock is attributed to the final, mixed permeate (mixed_permeate) exiting the module. The inlet and outlet StateBlocks of the permeate are used to only determine the permeate solute concentration for solvent and solute flux at the feed-side inlet and outlet, while the mixed_permeate StateBlock is used for mass balance based on the average flux.

Variables
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Water flux, solute flux", ":math:`j_v, j_s,_i`", "flux_mass_phase_comp", "[p]", ":math:`\text{kg/m}^2\text{/s}`"
   "Pore diffusivity of ion", ":math:`D_i,_p`", "diffus_pore_phase_comp", "[p]", ":math:`\text{m}^2\text{/s}`"
   "Convective hindrance factor", ":math:`k_i,_c`", "hindrance_factor_term_comp[[convective, diffusive],c]", "None", ":math:`\text{dimensionless}`"
   "Diffusive hindrance factor", ":math:`k_i,_d`", "hindrance_factor_term_comp[[convective, diffusive],c]", "None", ":math:`\text{dimensionless}`"
   "Pore Diffusivity of ion", ":math:`D_i,_p`", "diffus_pore_phase_comp", "[p]", ":math:`\text{m}^2\text{/s}`"
   "Pore radius", ":math:`r_p`", "radius_pore", "[p]", ":math:`\text{m}`"
   "Stokes radius", ":math:`r_s,_i`", "radius_stokes_comp[c]", "[p]", ":math:`\text{m}`"
   "rs/rp", ":math:`λ_i`", "lambda_comp[c]", "[p]", ":math:`\text{dimensionless}`"
   "Effective membrane thickness", ":math:`Δx_e`", "membrane_thickness_effective", "None", ":math:`\text{m}`"
   "Valency", ":math:`z_i`", "charge_comp[c]", "[p]", ":math:`\text{dimensionless}`"
   "Steric partitioning factor", ":math:`φ_i`", "partitioning_factor_steric_comp", "None", ":math:`\text{dimensionless}`"
   "Born solvation contribution to partitioning", ":math:`φ_{b_i}`", "partitioning_factor_born_solvation_comp", "[p]", ":math:`\text{dimensionless}`"
   "Gibbs free energy of solvation", ":math:`dG_{solv}`", "gibbs_solvation_comp", "None", ":math:`\text{J}`"
   "Membrane charge density", ":math:`c_x`", "membrane_charge_density", "None", ":math:`\text{mol/m}^3`"
   "Dielectric constant of medium (pore)", ":math:`Σ_p`", "dielectric_constant_pore", "[p]", ":math:`\text{dimensionless}`"
   "Dielectric constant of medium (feed) assumed equal to that of water", ":math:`Σ_f`", "dielectric_constant_feed", "[p]", ":math:`\text{dimensionless}`"
   "Concentration", ":math:`C_{i,j}`", "[feed,interface,pore_entrance,pore_exit,permeate].conc_mol)phase_comp", "[p,j]", ":math:`\text{kg/m}^3`"
   "Electric potential gradient between feed/interface", ":math:`xi`", "electric_potential_grad_feed_interface", "None", ":math:`\text{dimensionless}`"
   "Electronic charge", ":math:`e_o`", "electronic_charge", "[p]", ":math:`\text{C}`"
   "Absolute permittivity of vacuum", ":math:`Σ_o`", "vacuum_electric_permittivity", "None", ":math:`\text{F/m}`"
   "Boltzmann constant", ":math:`k_b`", "boltzmann_constant", "None", ":math:`\text{J/K}`"
   "Faraday's constant", ":math:`F`", "faraday_constant", "None", ":math:`\text{dimensionless}`"
   "Ideal gas constant", ":math:`R`", "gas_constant", "None", ":math:`\text{Check}`"

Relationships
---------------------------------------------------------------
.. csv-table::
   :header: "Description", "Equation"

   "Solvent flux in active layer (pore) domain", ":math:`J_s,_j = -D_{i,p}\frac{c_{i,j+1}-c_{i,j}}{δx_{j}}-0.5z_{i}(c_{i,j}+c_{i,j+1})D_{i,p}\frac{F}{RT}\frac{ψ_{j+1}-ψ_{j}}{δx_{j}}+0.5K_{i,c}(c_{i,j}+c_{i,j+1})J_{v}`"
   "Solute flux at feed/interface domain", ":math:`J_i = -k_{i}(C_{i,m}-C_{i,f})+J_{w}C_{i,m}-z_{i}C_{i,m}D_{i,∞}\frac{F}{RT}ξ`"
   "Solute flux - solvent flux relationship", ":math:`J_i = J_{v}c_{i,p}`" 
   "Diffusive hindered transport coefficient :math:`(λ_{i} ≤ 0.95)`", ":math:`K_{i,d} = \frac{1+(\frac{9}{8})λ_{i}ln(λ_{i})-1.56034λ_{i}+0.528155λ_{i}^{2}+1.91521λ_{i}^{3}-2.81903λ_{i}^{4}+0.270788λ_{i}^{5}-1.10115λ_{i}^{6}-0.435933λ_{i}^{7}}{(1-λ_{i})^{2}}`"
   "Diffusive hindered transport coefficient :math:`(λ_{i} > 0.95)`", ":math:`K_{i,d} = 0.984(\frac{1-λ_{i}}{λ_{i}})^{(5/2)}`"
   "Convective hindered transport coefficient", ":math:`K_{i,c} = \frac{1+3.867λ_{i}-1.907λ_{i}^{2}-0.834λ_{i}^{3}}{1+1.867λ_{i}-0.741λ_{i}^{2}}`"
   "Stokes pore radius ratio", ":math:`λ_{i} = \frac{r_{i,stokes}}{r_{pore}}`"
   "Pore diffusion coefficient", ":math:`D_{i,p} = K_{i,d}D_{i,∞}`"
   "Steric partitioning factor", ":math:`Φ_i = (1-λ_{i})^2`"
   "Born solvation partitioning", ":math:`Φ_b = exp(\frac{-ΔG_{i}}{k_{b}T})`"
   "Gibbs free energy of solvation", ":math:`ΔG = \frac{z_{i}^{2}e_{0}^{2}}{8πε_{0}r_{i}}(\frac{1}{ε_{pore}}-\frac{1}{ε_{f}})`"
   "Solvent flux (Hagen-Poiseuille)", ":math:`J_w = ΔP_{net}\frac{r_{pore}^{2}}{8vρ_{w}Δx_e} =((P_{f}-P_{p})-Δπ)\frac{r_{pore}^{2}}{8vρ_{w}Δx_e}`"
   "Membrane-solution interface equilibrium", ":math:`γ_{i,1}c_{i,1} = γ_{i,m}c_{i,m}Φ_{i}Φ_{b}exp(\frac{-z_{i}FΔψ_{D,m}}{RT})`"
   "Membrane-solution interface equilibrium", ":math:`γ_{i,N}c_{i,N} = γ_{i,p}c_{i,p}Φ_{i}Φ_{b}exp(\frac{-z_{i}FΔψ_{D,p}}{RT})`"
   

Scaling
-------
The DSPM-DE model includes support for scaling, such as providing default or calculating scaling factors for almost all variables.

   
Class Documentation
-------------------
* :mod:`watertap.unit_models.nanofiltration_DSPMDE_0D`


References
----------
Geraldes and Alves, 2008 `<https://doi.org/10.1016/j.memsci.2008.04.054>`_

Roy et al., 2015 `<http://dx.doi.org/10.1016/j.memsci.2015.06.030>`_

Labban et al., 2017 `<http://dx.doi.org/10.1016/j.memsci.2016.08.062>`_

Wang and Lin, 2021 `<https://doi.org/10.1016/j.memsci.2020.118809>`_
