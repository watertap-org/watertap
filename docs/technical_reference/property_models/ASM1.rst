ASM1 Property Package
=====================

This package implements property relationships for the treatment of wastewater using an activated sludge biological reactor as provided in `Henze, M. et al. (1987) <https://belinra.inrae.fr/doc_num.php?explnum_id=4467>`_.

This Activated Sludge Model no.1 (ASM1) property/reaction package:
   * supports 'H2O', 'S_I', 'S_S', 'X_I', 'X_S', 'X_BH', 'X_BA', 'X_P', 'S_O', 'S_NO', 'S_NH', 'S_ND', 'X_ND', and 'S_ALK' as components
   * supports only liquid phase

Sets
----
.. csv-table::
  :header: "Description", "Symbol", "Indices"

  "Components", ":math:`j`", "['H2O', 'S_I', 'S_S', 'X_I', 'X_S', 'X_BH', 'X_BA', 'X_P', 'S_O', 'S_NO', 'S_NH', 'S_ND', 'X_ND', 'S_ALK']"
  "Phases", ":math:`p`", "['Liq']"

Components
----------
.. csv-table::
  :header: "Description", "Symbol", "Variable"

  "Soluble inert organic matter, S_I", ":math:`S_I`", "S_I"
  "Readily biodegradable substrate S_S", ":math:`S_S`", "S_S"
  "Particulate inert organic matter, X_I", ":math:`X_I`", "X_I"
  "Slowly biodegradable substrate X_S", ":math:`X_S`", "X_S"
  "Active heterotrophic biomass X_B,H", ":math:`X_{B,H}`", "X_BH"
  "Active autotrophic biomass X_B,A", ":math:`X_{B,A}`", "X_BA"
  "Particulate products arising from biomass decay, X_P", ":math:`X_P`", "X_P"
  "Oxygen, S_O", ":math:`S_O`", "S_O"
  "Nitrate and nitrite nitrogen, S_NO", ":math:`S_{NO}`", "S_NO"
  "NH4 :math:`^{+}` + NH :math:`_{3}` Nitrogen, S_NH", ":math:`S_{NH}`", "S_NH"
  "Soluble biodegradable organic nitrogen, S_ND", ":math:`S_{ND}`", "S_ND"
  "Particulate biodegradable organic nitrogen, X_ND", ":math:`X_{ND}`", "X_ND"
  "Alkalinity, S_ALK", ":math:`S_{ALK}`", "S_ALK"

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Total volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"
   "Component mass concentrations", ":math:`C_j`", "conc_mass_comp", "[j]", ":math:`\text{kg/}\text{m}^3`"
   "Alkalinity in molar concentration", ":math:`A`", "alkalinity", "None", ":math:`\text{kmol HCO}_{3}^{-}\text{/m}^{3}`"

Stoichiometric Parameters
-------------------------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Value at 20 C", "Units"

   "Yield of cell COD formed per g N consumed, Y_A", ":math:`Y_A`", "Y_A", 0.24, ":math:`\text{dimensionless}`"
   "Yield of cell COD formed per g COD oxidized, Y_H", ":math:`Y_H`", "Y_H", 0.67, ":math:`\text{dimensionless}`"
   "Fraction of biomass yielding particulate products, f_p", ":math:`f_P`", "f_p", 0.08, ":math:`\text{dimensionless}`"
   "Mass fraction of N per COD in biomass, i_xb", ":math:`i_{XB}`", "i_xb", 0.08, ":math:`\text{dimensionless}`"
   "Mass fraction of N per COD in particulates, i_xp", ":math:`i_{XP}`", "i_xp", 0.06, ":math:`\text{dimensionless}`"

Kinetic Parameters
------------------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Value at 20 C", "Units"

   "Maximum specific growth rate for autotrophic biomass, mu_A", ":math:`µ_A`", "mu_A", 0.5, ":math:`\text{d}^{-1}`"
   "Maximum specific growth rate for heterotrophic biomass, mu_H", ":math:`µ_H`", "mu_H", 4.0, ":math:`\text{d}^{-1}`"
   "Half-saturation coefficient for heterotrophic biomass, K_S", ":math:`K_S`", "K_S", 0.01, ":math:`\text{kg COD/}\text{m}^{3}`"
   "Oxygen half-saturation coefficient for heterotrophic biomass, K_O,H", ":math:`K_{O,H}`", "K_OH", 0.0002, ":math:`\text{kg -COD/}\text{m}^{3}`"
   "Oxygen half-saturation coefficient for autotrophic biomass, K_O,A", ":math:`K_{O,A}`", "K_OA", 0.0004, ":math:`\text{kg -COD/}\text{m}^{3}`"
   "Nitrate half-saturation coefficient for denitrifying heterotrophic biomass, K_NO", ":math:`K_{NO}`", "K_NO", 0.0005, ":math:`\text{kg NO}_{3}\text{-N/}\text{m}^{3}`"
   "Decay coefficient for heterotrophic biomass, b_H", ":math:`b_H`", "b_H", 0.3, ":math:`\text{d}^{-1}`"
   "Decay coefficient for autotrophic biomass, b_A", ":math:`b_A`", "b_A", 0.05, ":math:`\text{d}^{-1}`"
   "Correction factor for mu_H under anoxic conditions, eta_g", ":math:`η_g`", "eta_g", 0.8, ":math:`\text{dimensionless}`"
   "Correction factor for hydrolysis under anoxic conditions, eta_h", ":math:`η_h`", "eta_h", 0.8, ":math:`\text{dimensionless}`"
   "Maximum specific hydrolysis rate, k_h", ":math:`k_h`", "k_h", 3.0, ":math:`\text{d}^{-1}`"
   "Half-saturation coefficient for hydrolysis of slowly biodegradable substrate, K_X", ":math:`K_X`", "K_X", 0.1, ":math:`\text{dimensionless}`"
   "Ammonia Half-saturation coefficient for autotrophic biomass, K_NH", ":math:`K_{NH}`", "K_NH", 1.0, ":math:`\text{kg NH}_{3}\text{-N/}\text{m}^{3}`"
   "Ammonification rate, k_a", ":math:`k_a`", "k_a", 0.00005, ":math:`\text{m}^{3}\text{/}\text{kg COD . d}`"

Properties
----------
.. csv-table::
  :header: "Description", "Symbol", "Variable", "Index", "Units"

  "Fluid specific heat capacity", ":math:`c_p`", "cp", "None", ":math:`\text{J/kg/K}`"
  "Mass density", ":math:`\rho`", "dens_mass", "[p]", ":math:`\text{kg/}\text{m}^3`"

Process Rate Equations
----------------------
.. csv-table::
   :header: "Description", "Equation"

   "Aerobic growth of heterotrophs", ":math:`ρ_1 = µ_{H}(\frac{S_{S}}{K_{S}+S_{S}})(\frac{S_{O}}{K_{O,H}+S_{O}})X_{B,H}`"
   "Anoxic growth of heterotrophs", ":math:`ρ_2 = µ_{H}(\frac{S_{S}}{K_{S}+S_{S}})(\frac{K_{O,H}}{K_{O,H}+S_{O}})(\frac{S_{NO}}{K_{NO}+S_{NO}})η_{g}X_{B,H}`"
   "Aerobic growth of autotrophs", ":math:`ρ_3 = µ_{A}(\frac{S_{NH}}{K_{NH}+S_{NH}})(\frac{S_{O}}{K_{O,A}+S_{O}})X_{B,A}`"
   "Decay of heterotrophs", ":math:`ρ_4 = b_{H}X_{B,H}`"
   "Decay of autotrophs", ":math:`ρ_5 = b_{H}X_{B,H}`"
   "Ammonification of soluble organic nitrogen", ":math:`ρ_6 = k_{a}S_{ND}X_{B,H}`"
   "Hydrolysis of entrapped organics", ":math:`ρ_7 = k_{H}(\frac{X_{S}/X_{B,H}}{K_{X}+(X_{S}/X_{B,H})})[(\frac{S_{O}}{K_{O,H}+S_{O}})+η_{h}(\frac{K_{O,H}}{K_{O,H}+S_{O}})(\frac{S_{NO}}{K_{NO}+S_{NO}})]X_{B,H}`"
   "Hydrolysis of entrapped organic nitrogen", ":math:`ρ_7 = k_{H}(\frac{X_{S}/X_{B,H}}{K_{X}+(X_{S}/X_{B,H})})[(\frac{S_{O}}{K_{O,H}+S_{O}})+η_{h}(\frac{K_{O,H}}{K_{O,H}+S_{O}})(\frac{S_{NO}}{K_{NO}+S_{NO}})]X_{B,H}(X_{ND}/X_{S})`"


Scaling
-------
Scaling for the ASM1 property package has yet to be implemented.

Class Documentation
-------------------
.. currentmodule:: watertap.property_models.activated_sludge.asm1_properties

.. autoclass:: ASM1ParameterBlock
    :members:
    :noindex:

.. autoclass:: ASM1ParameterData
    :members:
    :noindex:

.. autoclass:: _ASM1StateBlock
    :members:
    :noindex:

.. autoclass:: ASM1StateBlockData
    :members:
    :noindex:

.. currentmodule:: watertap.property_models.activated_sludge.asm1_reactions

.. autoclass:: ASM1ReactionParameterBlock
    :members:
    :noindex:

.. autoclass:: ASM1ReactionParameterData
    :members:
    :noindex:

.. autoclass:: _ASM1ReactionBlock
    :members:
    :noindex:

.. autoclass:: ASM1ReactionBlockData
    :members:
    :noindex:


References
----------
[1] Henze, M., Grady, C.P.L., Gujer, W., Marais, G.v.R., Matsuo, T.,
"Activated Sludge Model No. 1", 1987, IAWPRC Task Group on Mathematical Modeling
for Design and Operation of Biological Wastewater Treatment.
https://belinra.inrae.fr/doc_num.php?explnum_id=4467

[2] Alex, J. et al. Benchmark Simulation Model no.1 (BSM1). Lund University, 2008, 5-6.
https://www.iea.lth.se/publications/Reports/LTH-IEA-7229.pdf