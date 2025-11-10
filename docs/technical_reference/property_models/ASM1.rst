.. _ASM1:

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

  "Soluble inert organic matter", ":math:`S_I`", "S_I"
  "Readily biodegradable substrate", ":math:`S_S`", "S_S"
  "Particulate inert organic matter", ":math:`X_I`", "X_I"
  "Slowly biodegradable substrate", ":math:`X_S`", "X_S"
  "Active heterotrophic biomass", ":math:`X_{B,H}`", "X_BH"
  "Active autotrophic biomass", ":math:`X_{B,A}`", "X_BA"
  "Particulate products arising from biomass decay", ":math:`X_P`", "X_P"
  "Oxygen", ":math:`S_O`", "S_O"
  "Nitrate and nitrite nitrogen", ":math:`S_{NO}`", "S_NO"
  "NH4 :math:`^{+}` + NH :math:`_{3}` Nitrogen", ":math:`S_{NH}`", "S_NH"
  "Soluble biodegradable organic nitrogen", ":math:`S_{ND}`", "S_ND"
  "Particulate biodegradable organic nitrogen", ":math:`X_{ND}`", "X_ND"
  "Alkalinity", ":math:`S_{ALK}`", "S_ALK"

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

   "Yield of cell COD formed per g N consumed", ":math:`Y_A`", "Y_A", 0.24, ":math:`\text{dimensionless}`"
   "Yield of cell COD formed per g COD oxidized", ":math:`Y_H`", "Y_H", 0.67, ":math:`\text{dimensionless}`"
   "Fraction of biomass yielding particulate product", ":math:`f_p`", "f_p", 0.08, ":math:`\text{dimensionless}`"
   "Mass fraction of N per COD in biomass", ":math:`i_{XB}`", "i_xb", 0.08, ":math:`\text{dimensionless}`"
   "Mass fraction of N per COD in particulates", ":math:`i_{XP}`", "i_xp", 0.06, ":math:`\text{dimensionless}`"
   "Conversion factor applied for TSS calculation", ":math:`CODtoSS`", "CODtoSS", 0.75, ":math:`\text{dimensionless}`"
   "Conversion factor for BOD5 for raw wastewater", ":math:`BOD5_{factor, raw}`", "BOD5_factor_raw", 0.65, ":math:`\text{dimensionless}`"
   "Conversion factor for BOD5 for effluent", ":math:`BOD5_{factor, effluent}`", "BOD5_factor_effluent", 0.25, ":math:`\text{dimensionless}`"

Kinetic Parameters
------------------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Value at 20 C", "Units"

   "Maximum specific growth rate for autotrophic biomass", ":math:`µ_A`", "mu_A", 0.5, ":math:`\text{d}^{-1}`"
   "Maximum specific growth rate for heterotrophic biomass", ":math:`µ_H`", "mu_H", 4.0, ":math:`\text{d}^{-1}`"
   "Half-saturation coefficient for heterotrophic biomass", ":math:`K_S`", "K_S", 0.01, ":math:`\text{kg COD/}\text{m}^{3}`"
   "Oxygen half-saturation coefficient for heterotrophic biomass", ":math:`K_{O,H}`", "K_OH", 0.0002, ":math:`\text{kg -COD/}\text{m}^{3}`"
   "Oxygen half-saturation coefficient for autotrophic biomass", ":math:`K_{O,A}`", "K_OA", 0.0004, ":math:`\text{kg -COD/}\text{m}^{3}`"
   "Nitrate half-saturation coefficient for denitrifying heterotrophic biomass", ":math:`K_{NO}`", "K_NO", 0.0005, ":math:`\text{kg NO}_{3}\text{-N/}\text{m}^{3}`"
   "Decay coefficient for heterotrophic biomass", ":math:`b_H`", "b_H", 0.3, ":math:`\text{d}^{-1}`"
   "Decay coefficient for autotrophic biomass, ":math:`b_A`", "b_A", 0.05, ":math:`\text{d}^{-1}`"
   "Correction factor for mu_H under anoxic conditions", ":math:`η_g`", "eta_g", 0.8, ":math:`\text{dimensionless}`"
   "Correction factor for hydrolysis under anoxic conditions", ":math:`η_h`", "eta_h", 0.8, ":math:`\text{dimensionless}`"
   "Maximum specific hydrolysis rate", ":math:`k_h`", "k_h", 3.0, ":math:`\text{d}^{-1}`"
   "Half-saturation coefficient for hydrolysis of slowly biodegradable substrate", ":math:`K_X`", "K_X", 0.1, ":math:`\text{dimensionless}`"
   "Ammonia Half-saturation coefficient for autotrophic biomass", ":math:`K_{NH}`", "K_NH", 1.0, ":math:`\text{kg NH}_{3}\text{-N/}\text{m}^{3}`"
   "Ammonification rate", ":math:`k_a`", "k_a", 0.00005, ":math:`\text{m}^{3}\text{/}\text{kg COD . d}`"

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


Effluent Metrics
----------------
.. csv-table::
  :header: "Description", "Variable", "Default Regulatory Limit", "Equation"

  "Total suspended solids", ":math:`TSS`", ":math:`0.03 kg/m^{3}`", ":math:`TSS = CODtoSS(X_{S} + X_{I} + X_{BH} + X_{BA} + X_{P})`"
  "Five-day biological oxygen demand (raw wastewater)", ":math:`BOD5_{raw}`", ":math:`0.01 kg/m^{3}`", ":math:`BOD5_{raw} = BOD5_{factor, raw}(S_{S} + X_{S} + (1 - f_{p})(X_{BH} + X_{BA}))`"
  "Five-day biological oxygen demand (effluent)", ":math:`BOD5_{effluent}`", ":math:`0.01 kg/m^{3}`", ":math:`BOD5_{effluent} = BOD5_{factor, effluent}(S_{S} + X_{S} + (1 - f_{p})(X_{BH} + X_{BA}))`"
  "Chemical oxygen demand", ":math:`COD`", ":math:`0.1 kg/m^{3}`", ":math:`COD = S_{S} + S_{I} + X_{S} + X_{I} + X_{BH} + X_{BA} + X_{P}`"
  "Total Kjeldahl nitrogen", ":math:`TKN`", ":math:`None`", ":math:`TKN = S_{NH} + S_{ND} + X_{ND} + i_{XB}(X_{BH} + X_{BA}) + i_{XP}(X_{P} + X_{I})`"
  "Total nitrogen", ":math:`N_{total}`", ":math:`0.018 kg/m^{3}`", ":math:`N_{total} = TKN + S_{NO}`"

Class Documentation
-------------------
.. currentmodule:: watertap.property_models.unit_specific.activated_sludge.asm1_properties

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

.. currentmodule:: watertap.property_models.unit_specific.activated_sludge.asm1_reactions

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