ASM2D Property Package
======================

This package implements properties and reactions of an activated sludge model for biological nutrient removal from wastewater using an activated sludge biological reactor with biological phosphorus removal as provided in Henze, M. et al. (1999) [1].

This Activated Sludge Model no.2D (ASM2D) property/reaction package:
   * supports 'H2O', 'S_A', 'S_F', 'S_I', S_N2, S_NH4, S_NO3, S_O2, S_PO4, S_ALK, X_AUT, X_H, X_I, X_MeOH, X_MeP, X_PAO, X_PHA, X_PP, X_S, and X_TSS as components
   * supports only liquid phase

Limitations of the model, as noted in [1], are as follows:
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
  "Dinitrogen, N2. SN2 is assumed to be the only nitrogenous product of denitrification", ":math:`S_{N2}`", "S_N2"
  "Ammonium plus ammonia nitrogen.", ":math:`NH_4`", "NH_4"
  "Nitrate plus nitrite nitrogen (N03' + N02' -N). SN03 is assumed to include nitrate as well as nitrite nitrogen.", ":math:`S_{NO3}`", "S_NO3"
  "Dissolved oxygen", ":math:`S_{O2}`", "S_O2"
  "Inorganic soluble phosphorus, primarily ortho-phosphates.", ":math:`S_{PO4}`", "S_PO4"
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
   "Component mass concentrations", ":math:`C_j`", "conc_mass_comp", "[j]", ":math:`\text{kg/}\text{m}^3`"
   "Molar alkalinity", ":math:`A`", "alkalinity", "None", ":math:`\text{kmol HCO}_{3}^{-}\text{/m}^{3}`"

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
   "Saturation/inhibition coefficient for nitrate", ":math:`K_{NO3}`", "K_NO3", 0.0005, ":math:`\text{kg N/}\text{m}^{3}`"
   "Saturation coefficient for particulate COD", ":math:`K_{X}`", "K_X", 0.1, ":math:`\text{kg X_S/}\text{kg X_H}`"
   "Maximum growth rate on substrate", ":math:`µ_H`", "mu_H", 6, ":math:`\text{kg X_S/}\text{kg X_H . day}`"
   "Maximum rate for fermentation", ":math:`q_{fe}`", "q_fe", 3, ":math:`\text{kg S_F/}\text{kg X_H . day}`"
   "Rate constant for lysis and decay", ":math:`b_H`", "b_H", 0.4, ":math:`\text{day}^{-1}`"
   "Saturation coefficient for growth on SF", ":math:`K_F`", "K_F", 0.004, ":math:`\text{kg COD/}\text{m}^{3}`"
   "Saturation coefficient for fermentation of SF", ":math:`K_{fe}`", "K_fe", 0.004, ":math:`\text{d}^{-1}`"
   "Saturation coefficient for growth on acetate SA", ":math:`K_A`", "K_A", 0.004, ":math:`\text{kg COD/}\text{m}^{3}`"
   "Saturation coefficient for ammonium (nutrient)", ":math:`K_{NH4}`", "K_NH4", 0.00005, ":math:`\text{kg N/}\text{m}^{3}`"
   "Saturation coefficient for phosphate (nutrient)", ":math:`K_P`", "K_P", 0.00001, ":math:`\text{kg P/}\text{m}^{3}`"
   "Saturation coefficient for alkalinity (HCO3-)", ":math:`K_{ALK}`", "K_ALK", 0.0001, ":math:`\text{kmol HCO_{3}^{-}/}\text{m}^{3}`"
   "Rate constant for storage of X_PHA (base Xpp)", ":math:`q_{PHA}`", "q_PHA", 3, ":math:`\text{kg PHA/}\text{kg PAO . day}`"
   "Rate constant for storage of X_PP", ":math:`q_{PP}`", "q_PP", 1.5, ":math:`\text{kg PP/}\text{kg PAO . day}`"
   "Maximum growth rate of PAO", ":math:`µ_{PAO}`", "mu_PAO", 1, ":math:`\text{day}^{-1}`"
   "Rate for Lysis of X_PAO", ":math:`b_{PAO}`", "b_PAO", 0.2, ":math:`\text{day}^{-1}`"
   "Rate for Lysis of X_PP", ":math:`b_{PP}`", "b_PP", 0.2, ":math:`\text{day}^{-1}`"
   "Rate for Lysis of X_PHA", ":math:`b_{PHA}`", "b_PHA", 0.2, ":math:`\text{day}^{-1}`"
   "Saturation coefficient for phosphorus in storage of PP", ":math:`K_{PS}`", "K_PS", 0.0002, ":math:`\text{kg P/}\text{m}^3`"
   "Saturation coefficient for poly-phosphate", ":math:`K_{PP}`", "K_PP", 0.01, ":math:`\text{kg PP/}\text{kg PAO}`"
   "Maximum ratio of X_PP/X_PAO", ":math:`K_{MAX}`", "K_MAX", 0.34, ":math:`\text{kg PP/}\text{kg PAO}`"
   "Inhibition coefficient for PP storage", ":math:`K_{IPP}`", "K_IPP", 0.02, ":math:`\text{kg PP/}\text{kg PAO}`"
   "Saturation coefficient for PHA", ":math:`K_{PHA}`", "K_PHA", 0.01, ":math:`\text{kg PHA/}\text{kg PAO}`"
   "Maximum growth rate of X_AUT", ":math:`µ_{AUT}`", "mu_AUT", 1, ":math:`\text{day}^{-1}`"
   "Decay rate of X_AUT", ":math:`b_{AUT}`", "b_{AUT}", 0.15, ":math:`\text{day}^{-1}`"
   "Rate constant for P precipitation", ":math:`k_{PRE}`", "k_pre", 1000, ":math:`\text{m/}^{3}\text{kg Fe(OH)_3 . day}`"
   "Rate constant for redissolution", ":math:`k_{RED}`", "k_red", 0.6, ":math:`\text{day}^{-1}`"

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

   "Aerobic hydrolysis", ":math:`ρ_1 = K_{H}(\frac{S_{O2}}{K_{O2}+S_{O2}})(\frac{X_{S}/X_{H}}{K_{X}+X_{S}/X_{H}})X_{H}`"
   "Anoxic hydrolysis", ":math:`ρ_2 = K_{H}η_{NO3}(\frac{K_{O2}}{K_{O2}+S_{O2}})(\frac{S_{NO3}}{K_{NO3}+S_{NO3}})(\frac{X_{S}/X_{H}}{K_{X}+X_{S}/X_{H}})X_{H}`"
   "Anaerobic hydrolysis", ":math:`ρ_3 = K_{H}η_{fe}(\frac{K_{O2}}{K_{O2}+S_{O2}})(\frac{K_{NO3}}{K_{NO3}+S_{NO3}})(\frac{X_{S}/X_{H}}{K_{X}+X_{S}/X_{H}})X_{H}`"
   "Growth on fermentable substrates, S_F ", ":math:`ρ_4 = µ_{H}(\frac{S_{O2}}{K_{O2}+S_{O2}})(\frac{S_{F}}{K_{F}+S_{F}})(\frac{S_{F}}{S_{F}+S_{A}})(\frac{S_{NH4}}{K_{NH4}+S_{NH4}})(\frac{S_{PO4}}{K_{P}+S_{PO4}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})X_{H}`"
   "Growth on fermentation products, S_A", ":math:`ρ_5 = µ_{H}(\frac{S_{O2}}{K_{O2}+S_{O2}})(\frac{S_{A}}{K_{A}+S_{A}})(\frac{S_{A}}{S_{F}+S_{A}})(\frac{S_{NH4}}{K_{NH4}+S_{NH4}})(\frac{S_{PO4}}{K_{P}+S_{PO4}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})X_{H}`"
   "Denitrification with fermentable substrates, S_F", ":math:`ρ_6 = µ_{H}η_{NO3}(\frac{K_{O2}}{K_{O2}+S_{O2}})(\frac{S_{NO3}}{K_{NO3}+S_{NO3}})(\frac{S_{F}}{K_{F}+S_{F}})(\frac{S_{F}}{S_{F}+S_{A}})(\frac{S_{NH4}}{K_{NH4}+S_{NH4}})(\frac{S_{PO4}}{K_{P}+S_{PO4}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})X_{H}`"
   "Denitrification with fermentation products, S_A", ":math:`ρ_7 = µ_{H}η_{NO3}(\frac{K_{O2}}{K_{O2}+S_{O2}})(\frac{S_{NO3}}{K_{NO3}+S_{NO3}})(\frac{S_{A}}{K_{A}+S_{A}})(\frac{S_{A}}{S_{F}+S_{A}})(\frac{S_{NH4}}{K_{NH4}+S_{NH4}})(\frac{S_{PO4}}{K_{P}+S_{PO4}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})X_{H}`"
   "Fermentation", ":math:`ρ_8 = q_{fe}(\frac{K_{O2}}{K_{O2}+S_{O2}})(\frac{K_{NO3}}{K_{NO3}+S_{NO3}})(\frac{S_{F}}{K_{F}+S_{F}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})X_{H}`"
   "Lysis", ":math:`ρ_9 = b_{H}X_{H}`"
   "Storage of X_PHA", ":math:`ρ_{10} = q_{PHA}(\frac{S_{A}}{K_{A}+S_{A}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})(\frac{X_{PP}/X_{PAO}}{K_{PP}+X_{PP}/X_{PAO}})X_{PAO}`"
   "Aerobic storage of X_PP", ":math:`ρ_{11} = q_{PP}(\frac{S_{O2}}{K_{O2}+S_{O2}})(\frac{S_{PO4}}{K_{PS}+S_{PO4}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})(\frac{X_{PHA}/X_{PAO}}{K_{PHA}+X_{PHA}/X_{PAO}})(\frac{K_{MAX} - X_{PP}/X_{PAO}}{K_{IPP}+K_{MAX} - X_{PP}/X_{PAO}})X_{PAO}`"
   "Anoxic storage of X_PP", ":math:`ρ_{12} = ρ_{11}η_{NO3}(\frac{K_{O2}}{S_{O2}})(\frac{S_{NO3}}{K_{NO3}+S_{NO3}})`"
   "Aerobic growth on X_PHA", ":math:`ρ_{13} = µ_{PAO}(\frac{S_{O2}}{K_{O2}+S_{O2}})(\frac{S_{NH4}}{K_{NH4}+S_{NH4}})(\frac{S_{PO4}}{K_{P}+S_{PO4}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})(\frac{X_{PHA}/X_{PAO}}{K_{PHA}+X_{PHA}/X_{PAO}})X_{PAO}`"
   "Anoxic growth on X_PHA", ":math:`ρ_{14} = ρ_{13}η_{NO3}(\frac{K_{O2}}{S_{O2}})(\frac{S_{NO3}}{K_{NO3}+S_{NO3}})`"
   "Lysis of X_PAO", ":math:`ρ_{15} = b_{PAO}X_{PAO}(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})`"
   "Lysis of X_PP", ":math:`ρ_{16} = b_{PP}X_{PP}(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})`"
   "Lysis of X_PHA", ":math:`ρ_{17} = b_{PHA}X_{PHA}(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})`"
   "Aerobic growth of X_AUT", ":math:`ρ_{18} = µ_{AUT}(\frac{S_{O2}}{K_{O2}+S_{O2}})(\frac{S_{NH4}}{K_{NH4}+S_{NH4}})(\frac{S_{PO4}}{K_{P}+S_{PO4}})(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})X_{AUT}`"
   "Lysis of X_AUT", ":math:`ρ_{19} = k_{AUT}X_{AUT}`"
   "Precipitation of phosphorus with ferric hydroxide", ":math:`ρ_{20} = k_{PRE}S_{PO4}X_{MeOH}`"
   "Redissolution", ":math:`ρ_{21} = k_{RED}X_{MeP}(\frac{S_{ALK}}{K_{ALK}+S_{ALK}})`"


Scaling
-------
A thorough scaling routine for the ASM2D property package has yet to be implemented.

Class Documentation
-------------------
.. currentmodule:: watertap.property_models.activated_sludge.asm2d_properties

.. autoclass:: ASM2dParameterBlock
    :members:
    :noindex:

.. autoclass:: ASM2dParameterData
    :members:
    :noindex:

.. autoclass:: _ASM2dStateBlock
    :members:
    :noindex:

.. autoclass:: ASM2dStateBlockData
    :members:
    :noindex:

.. currentmodule:: watertap.property_models.activated_sludge.asm2d_reactions

.. autoclass:: ASM2dReactionParameterBlock
    :members:
    :noindex:

.. autoclass:: ASM2dReactionParameterData
    :members:
    :noindex:

.. autoclass:: _ASM2dReactionBlock
    :members:
    :noindex:

.. autoclass:: ASM2dReactionBlockData
    :members:
    :noindex:

References
----------
[1] M. Henze, W. Gujer, T. Mino, T. Matsuo, M.C. Wentzel, G. v. R. Marais, M.C.M. Van Loosdrecht, Activated sludge model No.2D, ASM2D, Water Science and Technology. 39 (1999) 165–182. https://doi.org/10.1016/S0273-1223(98)00829-4.


