ASM2D Property Package
======================

This package implements properties and reactions of an activated sludge model for biological nutrient removal from wastewater using an activated sludge biological reactor with biological phosphorus removal as provided in Henze, M. et al. (1999) [1].

This Activated Sludge Model no.2D (ASM2D) property/reaction package:
   * supports 'H2O', 'S_A', 'S_F', 'S_I', S_N2, S_NH4, S_NO3, S_O2, S_PO4, S_ALK, X_AUT, X_H, X_I, X_MeOH, X_MeP, X_PAO, X_PHA, X_PP, X_S, and X_TSS as components
   * supports only liquid phase

Limitations of the model, based on the original literature referenced, are as follows:
  * valid for municipal wastewater only
  * overflow of fermentation products (considered to be acetate) to the aeration tank cannot be modeled
  * influent wastewater requires sufficient concentrations of magnesium and potassium
  * pH should be near neutral
  * temperature in the range of 10-25°C

Sets
----
.. csv-table::
  :header: "Description", "Symbol", "Indices"

  "Components", ":math:`j`", "['H2O', 'S_A', 'S_F', 'S_I', S_N2, S_NH4, S_NO3, S_O2, S_PO4, S_ALK, X_AUT, X_H, X_I, X_MeOH, X_MeP, X_PAO, X_PHA, X_PP, X_S, X_TSS]"
  "Phases", ":math:`p`", "['Liq']"

Components
----------
The ASM2D model includes 19 components as outlined in the table below.

.. csv-table::
  :header: "Description", "Symbol", "Name in Model"

  "Fermentation products, considered to be acetate.", ":math:`S_A`", "S_A"
  "Fermentable, readily bio-degradable organic substrates", ":math:`S_F`", "S_F"
  "Inert soluble organic material.", ":math:`S_I`", "S_I"
  "Dinitrogen, N2. SN2 is assumed to be the only nitrogenous product of denitrification", ":math:`S_N_2`", "S_N2"
  "Ammonium plus ammonia nitrogen.", ":math:`NH_4`", "NH_4"
  "Nitrate plus nitrite nitrogen (N03' + N02' -N). SN03 is assumed to include nitrate as well as nitrite nitrogen.", ":math:`S_NO_3`", "S_NO3"
  "Dissolved oxygen", ":math:`S_O_2`", "S_O2"
  "Inorganic soluble phosphorus, primarily ortho-phosphates.", ":math:`S_PO4`", "S_PO4"
  "Alkalinity, [mol HCO_3- per m^3]", ":math:`S_{ALK}`", "S_ALK"
  "Autotrophic nitrifying organisms.", ":math:`X_{AUT}`", "X_AUT"
  "Heterotrophic organisms.", ":math:`X_H`", "X_H"
  "Inert particulate organic material.", ":math:`X_I`", "X_I"
  "Metal-hydroxides.", ":math:`X_{MeOH}`", "X_MeOH"
  "Metal-phosphate.", ":math:`X_{MeP}`", "X_MeP"
  "Phosphate-accumulating organisms.", ":math:`X_{PAO}`", "X_PAO"
  "A cell internal storage product of phosphorus-accumulating organisms, primarily comprising poly-hydroxy-alkanoates (PHA).", ":math:`X_{PHA}`", "X_PHA"
  "Poly-phosphate.", ":math:`X_{PP}`", "X_PP"
  "Slowly biodegradable substrates.", ":math:`X_S`", "X_S"
  "Total suspended solids, TSS.", ":math:`X_{TSS}`", "X_TSS"

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Total volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"
   "Component mass concentrations", ":math:`C_j`", "conc_mass_comp", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Molar alkalinity", ":math:`A`", "alkalinity", "[p]", ":math:`\text{kmol HCO}_{3}^{-}\text{/m}^{3}`"

Stoichiometric Coefficients
---------------------------
The table below provides typical values for stoichiometric coefficients (from Table 9 of reference).

.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Default Value", "Units"

   "N content of inert soluble COD S_I", ":math:`i_{NSI}`", "i_NSI", 0.01, ":math:`\text{dimensionless}`"
   "N content of fermentable substrate S_F", ":math:`i_{NSF}`", "i_NSF", 0.03, ":math:`\text{dimensionless}`"
   "N content of inert particulate COD X_I", ":math:`i_{NXI}`", "i_NXI", 0.02, ":math:`\text{dimensionless}`"
   "N content of slowly biodegradable substrate X_S", ":math:`i_{NXS}`", "i_NXS", 0.08, ":math:`\text{dimensionless}`"
   "N content of biomass, X_H, X_PAO, X_AUT", ":math:`i_{NBM}`", "i_NBM", 0.07, ":math:`\text{dimensionless}`"
   "P content of inert soluble COD S_I", ":math:`i_{PSI}`", "i_PSI", 0.00, ":math:`\text{dimensionless}`"
   "P content of fermentable substrate, S_F", ":math:`i_{SF}`", "i_SF", 0.01, ":math:`\text{dimensionless}`"
   "P content of inert particulate COD X_I", ":math:`i_{PXI}`", "i_PXI", 0.01, ":math:`\text{dimensionless}`"
   "P content of slowly biodegradable substrate X_S", ":math:`i_{PXS}`", "i_PXS", 0.01, ":math:`\text{dimensionless}`"
   "P content of biomass, X_H, X_PAO, X_AUT", ":math:`i_{PBM}`", "i_PBM", 0.02, ":math:`\text{dimensionless}`"
   "TSS to COD ratio for X_I", ":math:`i_{TSSXI}`", "i_TSSXI", 0.75, ":math:`\text{dimensionless}`"
   "TSS to COD ratio for X_S", ":math:`i_{TSSXS}`", "i_TSSXS", 0.75, ":math:`\text{dimensionless}`"
   "TSS to COD ratio for biomass, X_H, X_PAO, X_AUT", ":math:`i_{TSSBM}`", "i_TSSBM", 0.90, ":math:`\text{dimensionless}`"
   "Production of S_I in hydrolysis", ":math:`f_{SI}`", "f_SI", 0, ":math:`\text{dimensionless}`"
   "Yield coefficient for heterotrophic biomass X_H", ":math:`Y_{H}`", "Y_H", 0.625, ":math:`\text{dimensionless}`"
   "Fraction of inert COD generated in lysis", ":math:`f_{XI}`", "f_XI", 0.1, ":math:`\text{dimensionless}`"
   "Yield coefficient for P accumulating organisms (biomass/PHA)", ":math:`Y_{PAO}`", "Y_PAO", 0.625, ":math:`\text{dimensionless}`"
   "PP requirement (PO4 release) per PHA stored", ":math:`Y_{PO4}`", "Y_PO4", 0.40, ":math:`\text{dimensionless}`"
   "PHA requirement for PP storage", ":math:`Y_{PHA}`", "Y_PHA", 0, ":math:`\text{dimensionless}`"
   "Yield of autotrophic biomass per NO3- N", ":math:`Y_{A}`", "Y_A", 0.24, ":math:`\text{dimensionless}`"

Kinetic Parameters
------------------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Value at 20°C", "Units"

   "Hydrolysis rate constant", ":math:`K_H`", "K_H", 3, ":math:`\text{day}^{-1}`"
   "Anoxic hydrolysis reduction factor", ":math:`η_{NO3}`", "eta_NO3", 0.6, ":math:`\text{dimensionless}`"
   "Anaerobic hydrolysis reduction factor", ":math:`η_{fe}`", "eta_fe", 0.40, ":math:`\text{dimensionless}`"
   "Saturation/inhibition coefficient for oxygen", ":math:`K_{O2}`", "K_O2", 0.0002, ":math:`\text{kg O_2/}\text{m}^{3}`"
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

Relationships
-------------
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
A thorough scaling routine for the ASM2D property package has yet to be implemented.


References
----------
[1] M. Henze, W. Gujer, T. Mino, T. Matsuo, M.C. Wentzel, G. v. R. Marais, M.C.M. Van Loosdrecht, Activated sludge model No.2D, ASM2D, Water Science and Technology. 39 (1999) 165–182. https://doi.org/10.1016/S0273-1223(98)00829-4.


