.. _ASM2d_GHG:

ASM2D-GHG Property Package
===============================
.. raw:: html

    <style> .red {color:red} </style>
    <style> .lime {color:lime} </style>
    <style> .blue {color:blue} </style>

.. role:: red

.. role:: lime

.. role:: blue

This package is an extension of the :ref:`modified Activated Sludge Model No.2 (ASM2d) <modified_ASM2d>` and implements
properties and reactions of an activated sludge model for biological nutrient removal from wastewater using
an activated sludge biological reactor that tracks greenhouse gas emissions as provided
in Solis, B. et al. (2022) [1].

Throughout this documentation, text in :red:`red` has been removed in the ASM2d-N2O model, text in :lime:`lime` has
been added, and text in :blue:`blue` has been changed from the modified ASM2d implementation.

The following changes have been made to the modified ASM2d model as provided in [1]:
   * adds hydroxylamine (S_NH2OH), nitrous oxide (S_N2O), nitrogen monoxide (S_NO), and nitrite (S_NO2) as solutes
   * replaces X_AUT with X_AOB and X_NOB
   * adds additional process rate equations for the new solutes


This Activated Sludge Model no.2D  (ASM2D-GHG) property/reaction package with greenhouse gas extensions:
   * supports 'H2O', 'S_O2', 'S_F', 'S_A', S_I, S_NH4, S_HN2OH, S_N2O, S_NO, S_NO2, S_NO3, S_N2, S_PO4, S_IC, X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AOB, X_NOB, S_K, and S_Mg as components
   * supports only liquid phase


Sets
----
.. csv-table::
  :header: "Description", "Symbol", "Indices"

  "Components", ":math:`j`", "['H2O', 'S_O2', 'S_F', 'S_A', S_I, S_NH4, S_NH2OH, S_N2O, S_NO, S_NO2, S_NO3, S_N2, S_PO4, S_IC, X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AOB, X_NOB, S_K, S_Mg]"
  "Phases", ":math:`p`", "['Liq']"

Components
----------
The modified ASM2D model includes 18 components as outlined in the table below. :red:`Red` text indicates the component has been removed in the Modified ASM2d model, and :lime:`lime` text indicates the component has been added.


.. csv-table::
  :header: "Description", "Symbol", "Name in Model"

  "Dissolved oxygen", ":math:`S_{O2}`", "S_O2"
  "Fermentable, readily bio-degradable organic substrates", ":math:`S_F`", "S_F"
  "Fermentation products, considered to be acetate", ":math:`S_A`", "S_A"
  "Inert soluble organic material", ":math:`S_I`", "S_I"
  "Ammonium plus ammonia nitrogen", ":math:`S_{NH4}`", "S_NH4"
  ":lime:`Hydroxylamine`", ":math:`S_{NH2OH}`", "S_NH2OH"
  ":lime:`Nitrous oxide`", ":math:`S_{N2O}`", "S_N2O"
  ":lime:`Nitric oxide`", ":math:`S_{NO}`", "S_NO"
  ":lime:`Nitrite`", ":math:`S_{NO2}`", "S_NO2"
  ":blue:`Nitrate`", ":math:`S_{NO3}`", "S_NO3"
  "Dinitrogen", ":math:`S_{N2}`", "S_N2"
  "Inorganic soluble phosphorus, primarily ortho-phosphates", ":math:`S_{PO4}`", "S_PO4"
  "Inorganic carbon", ":math:`S_{IC}`", "S_IC"
  "Inert particulate organic material", ":math:`X_I`", "X_I"
  "Slowly biodegradable substrates", ":math:`X_S`", "X_S"
  "Heterotrophic organisms", ":math:`X_H`", "X_H"
  "Phosphate-accumulating organisms", ":math:`X_{PAO}`", "X_PAO"
  "Poly-phosphate", ":math:`X_{PP}`", "X_PP"
  "A cell internal storage product of phosphorus-accumulating organisms, primarily comprising poly-hydroxy-alkanoates (PHA)", ":math:`X_{PHA}`", "X_PHA"
  ":red:`Autotrophic nitrifying organisms`", ":math:`X_{AUT}`", "X_AUT"
  ":lime:`Ammonia oxidizing bacteria`", ":math:`X_{AOB}`", "X_AOB"
  ":lime:`Nitrite oxidizing bacteria`", ":math:`X_{NOB}`", "X_NOB"
  "Potassium", ":math:`S_{K}`", "S_K"
  "Magnesium", ":math:`S_{Mg}`", "S_Mg"

State variables
---------------
:red:`Red` text indicates the state variable has been removed in the Modified ASM2d model.

.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Total volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"
   "Component mass concentrations", ":math:`C_j`", "conc_mass_comp", "[j]", ":math:`\text{kg/}\text{m}^3`"

Stoichiometric Coefficients
---------------------------
:red:`Red` text indicates the stoichiometric coefficient has been removed in the Modified ASM2d model, :lime:`lime` text indicates the stoichiometric coefficient has been added, and :blue:`blue` text indicates the coefficient has had its value changed from its base ASM2d implementation.

:Note: Parameter values found in asm2dinit_bsm2.c: https://app.box.com/file/1812331425408

.. csv-table::
   :header: "Description", "Symbol", "Parameter", "Default Value", "Units"

   "C content of inert soluble COD S_I", ":math:`i_{CSI}`", "i_CSI", 0.36718, ":math:`\text{dimensionless}`"
   "C content of inert soluble COD S_F", ":math:`i_{CSF}`", "i_CSF", 0.31843, ":math:`\text{dimensionless}`"
   "C content of inert soluble COD S_A", ":math:`i_{CSA}`", "i_CSA", 0.375, ":math:`\text{dimensionless}`"
   "C content of inert soluble COD X_I", ":math:`i_{CXI}`", "i_CXI", 0.36178, ":math:`\text{dimensionless}`"
   "C content of inert soluble COD X_S", ":math:`i_{CXS}`", "i_CXS", 0.31843, ":math:`\text{dimensionless}`"
   ":blue:`C content of inert soluble COD X_B`", ":math:`i_{CBM}`", "i_CBM", 0.36612, ":math:`\text{dimensionless}`"
   ":lime:`C content of inert soluble COD X_PHA`", ":math:`i_{CXPHA}`", "i_CXPHA", 0.3, ":math:`\text{dimensionless}`"
   "N content of inert soluble COD S_I", ":math:`i_{NSI}`", "i_NSI", 0.06003, ":math:`\text{dimensionless}`"
   "N content of fermentable substrate S_F", ":math:`i_{NSF}`", "i_NSF", 0.03552, ":math:`\text{dimensionless}`"
   "N content of inert particulate COD X_I", ":math:`i_{NXI}`", "i_NXI", 0.06003, ":math:`\text{dimensionless}`"
   "N content of slowly biodegradable substrate X_S", ":math:`i_{NXS}`", "i_NXS", 0.03552, ":math:`\text{dimensionless}`"
   "N content of biomass, X_H, X_PAO, X_AOB, X_NOB", ":math:`i_{NBM}`", "i_NBM", 0.08615, ":math:`\text{dimensionless}`"
   ":lime:`P content of inert soluble COD S_I`", ":math:`i_{PSI}`", "i_PSI", 0.00649, ":math:`\text{dimensionless}`"
   "P content of fermentable substrate, S_F", ":math:`i_{PSF}`", "i_PSF", 0.00559, ":math:`\text{dimensionless}`"
   "P content of inert particulate COD X_I", ":math:`i_{PXI}`", "i_PXI", 0.00649, ":math:`\text{dimensionless}`"
   "P content of slowly biodegradable substrate X_S", ":math:`i_{PXS}`", "i_PXS", 0.00559, ":math:`\text{dimensionless}`"
   "P content of biomass, X_H, X_PAO, X_AOB, X_NOB", ":math:`i_{PBM}`", "i_PBM", 0.02154, ":math:`\text{dimensionless}`"
   "TSS to COD ratio for X_I", ":math:`i_{TSSXI}`", "i_TSSXI", 0.75, ":math:`\text{dimensionless}`"
   "TSS to COD ratio for X_S", ":math:`i_{TSSXS}`", "i_TSSXS", 0.75, ":math:`\text{dimensionless}`"
   "TSS to COD ratio for biomass, X_H, X_PAO, X_AUT", ":math:`i_{TSSBM}`", "i_TSSBM", 0.90, ":math:`\text{dimensionless}`"
   "Production of S_I in hydrolysis", ":math:`f_{SI}`", "f_SI", 0, ":math:`\text{dimensionless}`"
   "Yield coefficient for heterotrophic biomass X_H", ":math:`Y_{H}`", "Y_H", 0.625, ":math:`\text{dimensionless}`"
   ":red:`Fraction of inert COD generated in lysis`", ":math:`f_{XI}`", "f_XI", 0.1, ":math:`\text{dimensionless}`"
   "Yield coefficient for P accumulating organisms (biomass/PHA)", ":math:`Y_{PAO}`", "Y_PAO", 0.625, ":math:`\text{dimensionless}`"
   ":blue:`PP requirement (PO4 release) per PHA stored`", ":math:`Y_{PO4}`", "Y_PO4", 0.4, ":math:`\text{dimensionless}`"
   "PHA requirement for PP storage", ":math:`Y_{PHA}`", "Y_PHA", 0.2, ":math:`\text{dimensionless}`"
   ":red:`Yield of autotrophic biomass per NO3- N`", ":math:`Y_{A}`", "Y_A", 0.24, ":math:`\text{dimensionless}`"
   ":lime:`Yield of ammonia oxidizing bacteria`", ":math:`Y_{AOB}`", "Y_AOB", 0.18, ":math:`\text{dimensionless}`"
   ":lime:`Yield of nitrite oxidizing bacteria`", ":math:`Y_{NOB}`", "Y_NOB", 0.08, ":math:`\text{dimensionless}`"
   "Potassium coefficient for polyphosphates", ":math:`i_{KXPP}`", "i_KXPP", 0.4204, ":math:`\text{dimensionless}`"
   "Magnesium coefficient for polyphosphates", ":math:`i_{MgXPP}`", "i_MgXPP", 0.2614, ":math:`\text{dimensionless}`"
   ":lime:`Anoxic growth factor`", ":math:`nG`", "nG", 1, ":math:`\text{dimensionless}`"

Kinetic Parameters
------------------
:red:`Red` text indicates the parameter has been removed in the Modified ASM2d model, :lime:`lime` text indicates the parameter has been added, and :blue:`blue` text indicates the parameter has had its value changed from its base ASM2d implementation.

.. csv-table::
   :header: "Description", "Symbol", "Parameter", "Value at 20°C", "Units"

   "Hydrolysis rate constant", ":math:`K_H`", "K_H", 2.46, ":math:`\text{day}^{-1}`"
   "Anoxic hydrolysis reduction factor for nitrate", ":math:`hl_{NO3}`", "hl_NO3", 0.6, ":math:`\text{dimensionless}`"
   "Anaerobic hydrolysis reduction factor", ":math:`hl_{fe}`", "hl_fe", 0.40, ":math:`\text{dimensionless}`"
   "Saturation/inhibition coefficient for oxygen", ":math:`KL_{O2}`", "KL_O2", 0.0002, ":math:`\text{kg O_2/}\text{m}^{3}`"
   "Saturation/inhibition coefficient for nitrate", ":math:`KL_{NO3}`", "KL_NO3", 0.0005, ":math:`\text{kg N/}\text{m}^{3}`"
   "Saturation coefficient for particulate COD", ":math:`KL_{X}`", "KL_X", 0.1, ":math:`\text{kg X_S/}\text{kg X_H}`"
   ":lime:`Anoxic hydrolysis reduction factor for nitrite`", ":math:`hl_{NO2}`", "hl_NO2", 0.6, ":math:`\text{dimensionless}`"
   ":lime:`Saturation/inhibition coefficient for nitrite`", ":math:`KL_{NO2}`", "KL_NO2", 0.0005, ":math:`\text{kg N/}\text{m}^{3}`"
   "Maximum growth rate on substrate", ":math:`µ_H`", "mu_H", 4.23, ":math:`\text{kg X_S/}\text{kg X_H . day}`"
   "Maximum rate for fermentation", ":math:`q_{fe}`", "q_fe", 2.11, ":math:`\text{kg S_F/}\text{kg X_H . day}`"
   ":blue:`Denitrification reduction factor for nitrate to nitrite`", ":math:`hH_{NO3}`", "hH_NO3", 0.28, ":math:`\text{dimensionless}`"
   "Rate constant for lysis and decay", ":math:`b_H`", "b_H", 0.28, ":math:`\text{day}^{-1}`"
   "Anoxic reduction factor for denitrification", ":math:`hH_{NO3, end}`", "hH_NO3_end", 0.5, ":math:`\text{dimensionless}`"
   "Saturation coefficient for growth on SF", ":math:`K_F`", "K_F", 0.004, ":math:`\text{kg COD/}\text{m}^{3}`"
   "Saturation coefficient for fermentation of SF", ":math:`K_{fe}`", "K_fe", 0.004, ":math:`\text{d}^{-1}`"
   "Saturation coefficient for growth on acetate SA", ":math:`KH_A`", "KH_A", 0.004, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":blue:`Saturation/inhibition coefficient for oxygen`", ":math:`KH_{O2}`", "KH_O2", 0.0001, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation/inhibition coefficient for nitrate`", ":math:`KH_{NO3}`", "KH_NO3", 0.0005, ":math:`\text{kg COD/}\text{m}^{3}`"
   "Saturation coefficient for ammonium (nutrient)", ":math:`KH_{NH4}`", "KH_NH4", 0.00005, ":math:`\text{kg N/}\text{m}^{3}`"
   "Saturation coefficient for phosphate (nutrient)", ":math:`KH_{PO4}`", "KH_PO4", 0.00001, ":math:`\text{kg P/}\text{m}^{3}`"
   ":lime:`Denitrification reduction factor for nitrite to nitric oxide`", ":math:`hH_{NO2}`", "hH_NO2", 0.16, ":math:`\text{dimensionless}`"
   ":lime:`Denitrification reduction factor for nitric oxide to nitrous oxide`", ":math:`hH_{NO}`", "hH_NO", 0.35, ":math:`\text{dimensionless}`"
   ":lime:`Denitrification reduction factor for nitrous oxide to dinitrogen`", ":math:`hH_{N2O}`", "hH_N2O", 0.35, ":math:`\text{dimensionless}`"
   ":lime:`Saturation/inhibition coefficient for oxygen (NO3 to NO2)`", ":math:`KH2_{O2}`", "KH2_O2", 0.0001, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation/inhibition coefficient for oxygen (NO2 to NO)`", ":math:`KH3_{O2}`", "KH3_O2", 0.0001, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation/inhibition coefficient for oxygen (NO to N2O)`", ":math:`KH4_{O2}`", "KH4_O2", 0.0001, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation/inhibition coefficient for oxygen (N2O to N2)`", ":math:`KH5_{O2}`", "KH5_O2", 0.0001, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation coefficient for growth on SF (NO3 to NO2)`", ":math:`K_{F2}`", "K_F2", 0.02, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation coefficient for growth on SF (NO2 to NO)`", ":math:`K_{F3}`", "K_F3", 0.02, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation coefficient for growth on SF (NO to N2O)`", ":math:`K_{F4}`", "K_F4", 0.02, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation coefficient for growth on SF (N2O to N2)`", ":math:`K_{F5}`", "K_F5", 0.04, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation coefficient for growth on acetate SA (NO3 to NO2)`", ":math:`KH_{A2}`", "KH_A2", 0.02, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation coefficient for growth on acetate SA (NO2 to NO)`", ":math:`KH_{A3}`", "KH_A3", 0.02, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation coefficient for growth on acetate SA (NO to N2O)`", ":math:`KH_{A4}`", "KH_A4", 0.02, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation coefficient for growth on acetate SA (N2O to N2)`", ":math:`KH_{A5}`", "KH_A5", 0.04, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation/inhibition coefficient for nitrite`", ":math:`KH_{NO2}`", "KH_NO2", 0.0002, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation/inhibition coefficient for nitric oxide`", ":math:`KH_{NO}`", "KH_NO", 0.00005, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Saturation/inhibition coefficient for nitrous oxide`", ":math:`KH_{N2O}`", "KH_N2O", 0.00005, ":math:`\text{kg COD/}\text{m}^{3}`"
   ":lime:`Inhibition coefficient for nitric oxide`", ":math:`KI_{NO}`", "KI_NO", 0.0003, ":math:`\text{kg COD/}\text{m}^{3}`"
   "Rate constant for storage of X_PHA", ":math:`q_{PHA}`", "q_PHA", 2.46, ":math:`\text{kg PHA/}\text{kg PAO . day}`"
   "Rate constant for storage of X_PP", ":math:`q_{PP}`", "q_PP", 1.23, ":math:`\text{kg PP/}\text{kg PAO . day}`"
   "Maximum growth rate of PAO", ":math:`µ_{PAO}`", "mu_PAO", 0.82, ":math:`\text{day}^{-1}`"
   "Rate for Lysis of X_PAO", ":math:`b_{PAO}`", "b_PAO", 0.14, ":math:`\text{day}^{-1}`"
   "Rate for Lysis of X_PP", ":math:`b_{PP}`", "b_PP", 0.14, ":math:`\text{day}^{-1}`"
   "Rate for Lysis of X_PHA", ":math:`b_{PHA}`", "b_PHA", 0.14, ":math:`\text{day}^{-1}`"
   "Saturation coefficient for oxygen", ":math:`KP_{O2}`", "KP_O2", 0.0002, ":math:`\text{kg P/}\text{m}^3`"
   "Saturation coefficient for nitrate", ":math:`KP_{NO3}`", "KP_NO3", 0.0005, ":math:`\text{kg P/}\text{m}^3`"
   "Saturation coefficient for acetate", ":math:`KP_A`", "KP_A", 0.004, ":math:`\text{kg P/}\text{m}^3`"
   "Saturation coefficient for ammonium", ":math:`KP_{NH4}`", "KP_NH4", 0.00005, ":math:`\text{kg P/}\text{m}^3`"
   "Saturation coefficient for phosphate for growth", ":math:`KP_{PO4}`", "KP_PO4", 0.00001, ":math:`\text{kg P/}\text{m}^3`"
   "Saturation coefficient for phosphate for XPP formation", ":math:`KP_P`", "KP_P", 0.0002, ":math:`\text{kg P/}\text{m}^3`"
   "Saturation coefficient for poly-phosphate", ":math:`KP_{PP}`", "KP_PP", 0.01, ":math:`\text{kg PP/}\text{kg PAO}`"
   "Maximum ratio of X_PP/X_PAO", ":math:`K_{MAX}`", "K_MAX", 0.34, ":math:`\text{kg PP/}\text{kg PAO}`"
   "Inhibition coefficient for PP storage", ":math:`KI_{PP}`", "KI_PP", 0.02, ":math:`\text{kg PP/}\text{kg PAO}`"
   "Saturation coefficient for PHA", ":math:`KP_{PHA}`", "KP_PHA", 0.01, ":math:`\text{kg PHA/}\text{kg PAO}`"
   ":red:`Maximum growth rate of X_AUT`", ":math:`µ_{AUT}`", "mu_AUT", 0.61, ":math:`\text{day}^{-1}`"
   ":red:`Decay rate of X_AUT`", ":math:`b_{AUT}`", "b_AUT", 0.09, ":math:`\text{day}^{-1}`"
   ":blue:`Reduction factor for denitrification`", ":math:`hH_{NO3}`", "hH_NO3", 0.28, ":math:`\text{dimensionless}`"
   "Anoxic reduction factor for endogenous respiration", ":math:`hH_{NO3, end}`", "hH_NO3_end", 0.5, ":math:`\text{dimensionless}`"
   ":blue:`Reduction factor under anoxic conditions`", ":math:`hP_{NO3}`", "hP_NO3", 0.28, ":math:`\text{dimensionless}`"
   "Anoxic reduction factor for decay of PAOs", ":math:`hP_{NO3, end}`", "hP_NO3_end", 0.33, ":math:`\text{dimensionless}`"
   "Anoxic reduction factor for decay of PP", ":math:`hPP_{NO3, end}`", "hPP_NO3_end", 0.33, ":math:`\text{dimensionless}`"
   "Anoxic reduction factor for decay of PHA", ":math:`hPHA_{NO3, end}`", "hPHA_NO3_end", 0.33, ":math:`\text{dimensionless}`"
   ":red:`Anoxic reduction factor for decay of autotrophs`", ":math:`hAUT_{NO3, end}`", "hAUT_NO3_end", 0.33, ":math:`\text{dimensionless}`"
   ":lime:`Reduction factor for denitrification (nitrite to nitric oxide)`", ":math:`hP_{NO2}`", "hP_NO2", 0.16, ":math:`\text{dimensionless}`"
   ":lime:`Reduction factor for denitrification (nitric oxide to nitrous oxide)`", ":math:`hP_{NO}`", "hP_NO", 0.35, ":math:`\text{dimensionless}`"
   ":lime:`Reduction factor for denitrification (nitrous oxide to dinitrogen)`", ":math:`hP_{N2O}`", "hP_N2O", 0.35, ":math:`\text{dimensionless}`"
   ":lime:`Saturation/inhibition coefficient for oxygen (NO3 to NO2)`", ":math:`KP2_{O2}`", "KP2_O2", 0.0001, ":math:`\text{kg P/}\text{m}^3`"
   ":lime:`Saturation/inhibition coefficient for oxygen (NO2 to NO)`", ":math:`KP3_{O2}`", "KP3_O2", 0.0001, ":math:`\text{kg P/}\text{m}^3`"
   ":lime:`Saturation/inhibition coefficient for oxygen (NO to N2O)`", ":math:`KP4_{O2}`", "KP4_O2", 0.0001, ":math:`\text{kg P/}\text{m}^3`"
   ":lime:`Saturation/inhibition coefficient for oxygen (N20 to N2)`", ":math:`KP5_{O2}`", "KP5_O2", 0.0001, ":math:`\text{kg P/}\text{m}^3`"
   ":lime:`Saturation/inhibition coefficient for nitrite`", ":math:`KP_{NO2}`", "KP_NO2", 0.0002, ":math:`\text{kg P/}\text{m}^3`"
   ":lime:`Saturation/inhibition coefficient for nitric oxide`", ":math:`KP_{NO}`", "KP_NO", 0.00005, ":math:`\text{kg P/}\text{m}^3`"
   ":lime:`Saturation/inhibition coefficient for nitrous oxide`", ":math:`KP_{N2O}`", "KP_N2O", 0.00005, ":math:`\text{kg P/}\text{m}^3`"
   :NOTE: Need to double check the units on the variables below and remove if unused
   ":lime:`Maximum rate for the AMO reaction`", ":math:`q_{AOB, AMO}`", "q_AOB,AMO", 5.2, ":math:`\text{kg AOB/}\text{kg PAO . day}`"
   ":lime:`Maximum AOB growth rate`", ":math:`µ_{AOB, HAO}`", "mu_AOB,HAO", 0.61, ":math:`\text{day}^{-1}`"
   ":lime:`Maximum rate for the HAO reaction`", ":math:`q_{AOB, HAO}`", "q_AOB,HAO", 5.2, ":math:`\text{kg AOB/}\text{kg PAO . day}`"
   ":lime:`Maximum N2O production rate by NH2OH oxidation pathway`", ":math:`q_{AOB, N2O, NN}`", "q_AOB,N2O,NN", 0.0078, ":math:`\text{kg AOB/}\text{kg PAO . day}`"
   ":lime:`Maximum N2O production rate by ND pathway`", ":math:`q_{AOB, N2O, ND}`", "q_AOB,N2O,ND", 1.3008, ":math:`\text{kg AOB/}\text{kg PAO . day}`"
   ":lime:`Decay rate of AOB`", ":math:`b_{AOB}`", "b_AOB", 0.096, ":math:`\text{day}^{-1}`"
   ":lime:`Anoxic reduction factor for AOB decay`", ":math:`h_{AOB, NO3, end}`", "h_AOB,NO3,end", 0.33, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for oxygen (AMO reaction)`", ":math:`KAOB1_{O2}`", "KAOB1_O2", 1, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for NH4`", ":math:`KAOB_{NH4}`", "KAOB_NH4", 0.004, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for oxygen (AMO reaction)`", ":math:`KAOB_{O2}`", "KAOB2_O2", 0.6, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for NH2OH`", ":math:`KAOB_{NH2OH}`", "KAOB_NH2OH", 0.3, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for phosphate`", ":math:`KAOB_{P}`", "KAOB_P", 0.01, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for NO (from HAO)`", ":math:`KAOB_{HAO, NO}`", "KAOB_HAO,NO", 0.0003, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for NO (from NirK)`", ":math:`KAOB_{NN, NO}`", "KAOB_NN,NO", 0.008, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for free nitrous acid`", ":math:`KAOB_{HNO2}`", "KAOB_HNO2", 0.0006, ":math:`\text{dimensionless}`"
   ":lime:`AOB constant for O2 effect on the ND pathway`", ":math:`KAOB_{ND, O2}`", "KAOB_ND,O2", 0.5, ":math:`\text{dimensionless}`"
   ":lime:`N2O constant for production inhibition by O2`", ":math:`KAOB_{I, O2}`", "KAOB_I,O2", 0.8, ":math:`\text{dimensionless}`"
   ":lime:`AOB affinity constant for nitrate`", ":math:`KAOB_{NO3}`", "KAOB_NO3", 0.5, ":math:`\text{dimensionless}`"
   ":lime:`Maximum NOB growth rate`", ":math:`µ_{NOB}`", "mu_NOB", 0.61, ":math:`\text{day}^{-1}`"
   ":lime:`Decay rate of  NOB`", ":math:`b_{NOB}`", "b_NOB", 0.096, ":math:`\text{day}^{-1}`"
   ":lime:`Anoxic reduction factor for NOB decay`", ":math:`h_{NOB, NO3, end}`", "h_NOB,NO3,end", 0.33, ":math:`\text{dimensionless}`"
   ":lime:`NOB affinity constant for oxygen`", ":math:`K_{NOB, O2}`", "k_NOB,O2", 1.2, ":math:`\text{dimensionless}`"
   ":lime:`NOB affinity constant for nitrite`", ":math:`K_{NOB, NO2}`", "k_NOB,NO2", 1e-6, ":math:`\text{dimensionless}`"
   ":lime:`NOB affinity constant for phosphate`", ":math:`K_{NOB, P}`", "k_NOB,P", 0.01, ":math:`\text{dimensionless}`"
   ":lime:`NOB affinity constant for nitrate`", ":math:`K_{NOB, NO3}`", "k_NOB,NO3", 0.5, ":math:`\text{dimensionless}`"
    :NOTE: Add additional parameters if anything else is used

Properties
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Fluid specific heat capacity", ":math:`c_p`", "cp", "None", ":math:`\text{J/kg/K}`"
   "Mass density", ":math:`\rho`", "dens_mass", "[p]", ":math:`\text{kg/}\text{m}^3`"

Process Rate Equations
----------------------
Equations marked "(with decay)" indicate that the decay of heterotrophs and autotrophs is dependent on the electron acceptor present. Equations marked "(without decay)" indicate that the decay of heterotrophs and autotrophs does not change.

:Note: Process rate equations found in asm2d_PSFe_GHG.c: https://app.box.com/file/1812340348414

.. csv-table::
   :header: "Description", "Expression"

   "MonodX", ":math:`(\frac{X_{S}/X_{H}}{KL_{X}+X_{S}/X_{H}})`"
   "MonodPHA", ":math:`(\frac{X_{PHA}/(X_{PAO}+1e-6)}{KP_{PHA}+X_{PHA}/(X_{PAO}+1e-6)})`"
   "MonodIPP", ":math:`(\frac{K_{MAX} - X_{PP}/(X_{PAO}+1e-6)}{K_{IPP}+K_{MAX} - X_{PP}/(X_{PAO}+1e-6)})`"
   "MonodPP", ":math:`(\frac{X_{PP}/(X_{PAO}+1e-6)}{KP_{PP}+X_{PP}/(X_{PAO}+1e-6)})`"
   "MonodPO4", ":math:`(\frac{S_{PO4}}{KH_{PO4}+S_{PO4}})`"
   "MonodNH4", ":math:`(\frac{S_{NH4}}{KH_{NH4}+S_{NH4}})`"
   "MonodSPO4_P", ":math:`(\frac{S_{PO4}}{KP_{P}+S_{PO4}})`"
   "MonodSNH4_P", ":math:`(\frac{S_{NH4}}{KP_{NH4}+S_{NH4}})`"
   "fSO2", ":math:`(/frac{S_{O2}}{KAOB_{ND,O2}+(1-2(/frac{KAOB_{ND,O2}}{KAOB_{I_O2}})^(0.5)S_{O2})+(/frac{S_{O2}^{2}}{KAOB_{I,O2}})})`"

.. csv-table::
   :header: "Description", "Equation"

   "Process 1", ":math:`ρ_1 = K_{H}(\frac{S_{O2}}{KL_{O2}+S_{O2}})MonodX*X_{H}`"
   "Process 2", ":math:`ρ_2 = K_{H}hl_{NO3}(\frac{KL_{O2}}{KL_{O2}+S_{O2}S_{NO3}/(KL_{NO3}+S_{NO3})})MonodX*X_{H}`"
   "Process 3", ":math:`ρ_3 = K_{H}hl_{NO2}(\frac{KL_{O2}}{KL_{O2}+S_{O2}S_{NO2}/(KL_{NO3}+S_{NO2})})MonodX*X_{H}`"
   "Process 4", ":math:`ρ_4 = K_{H}hl_{fe}(\frac{KL_{O2}}{KL_{O2}+S_{O2}})(/frac{KL_{NO3}}{KL_{NO3}+S_{NO3}+S_{NO2}})MonodX*X_{H}`"
   "Process 5", ":math:`ρ_5 = µ_{H}(\frac{S_{O2}}{KH_{O2}+S_{O2}})(\frac{S_{F}}{K_{F}+S_{F}})(\frac{S_{F}}{S_{F}+S_{A}})MonodNH4*MonodPO4*X_{H}`"
   "Process 6", ":math:`ρ_6 = µ_{H}(\frac{S_{O2}}{KH_{O2}+S_{O2}})(\frac{S_{A}}{KH_{A}+S_{A}})(\frac{S_{A}}{S_{F}+S_{A}})MonodNH4*MonodPO4*X_{H}`"
   "Process 7", ":math:`ρ_7 = µ_{H}hH_{NO3}(\frac{KH2_{O2}}{KH2_{O2}+S_{O2}})(\frac{S_{F}}{K_{F2}+S_{F}})(\frac{S_{F}}{S_{F}+S_{A}})(\frac{S_{NO3}}{KH_{NO3}+S_{NO3}})MonodNH4*MonodPO4*X_{H}`"
   "Process 8", ":math:`ρ_8 = µ_{H}hH_{NO2}(\frac{KH3_{O2}}{KH3_{O2}+S_{O2}})(\frac{S_{F}}{K_{F3}+S_{F}})(\frac{S_{F}}{S_{F}+S_{A}})(\frac{S_{NO2}}{KH_{NO2}+S_{NO2}})MonodNH4*MonodPO4*X_{H}`"
   "Process 9", ":math:`ρ_9 = µ_{H}hH_{NO}(\frac{KH4_{O2}}{KH4_{O2}+S_{O2}})(\frac{S_{F}}{K_{F4}+S_{F}})(\frac{S_{F}}{S_{F}+S_{A}})(\frac{S_{NO}}{KH_{NO}+S_{NO}+S_{NO}^{2}/KI_{NO}})MonodNH4*MonodPO4*X_{H}`"
   "Process 10", ":math:`ρ_10 = µ_{H}hH_{N2O}(\frac{KH5_{O2}}{KH5_{O2}+S_{O2}})(\frac{S_{F}}{K_{F5}+S_{F}})(\frac{S_{F}}{S_{F}+S_{A}})(\frac{S_{N2O}}{KH_{N2O}+S_{N2O}})MonodNH4*MonodPO4*X_{H}`"
   "Process 11", ":math:`ρ_11 = µ_{H}hH_{NO3}(\frac{KH2_{O2}}{KH2_{O2}+S_{O2}})(\frac{S_{A}}{KH_{A2}+S_{A}})(\frac{S_{A}}{S_{F}+S_{A}})(\frac{S_{NO3}}{KH_{NO3}+S_{NO3}})MonodNH4*MonodPO4*X_{H}`"
   "Process 12", ":math:`ρ_12 = µ_{H}hH_{NO2}(\frac{KH3_{O2}}{KH3_{O2}+S_{O2}})(\frac{S_{A}}{KH_{A3}+S_{A}})(\frac{S_{A}}{S_{F}+S_{A}})(\frac{S_{NO2}}{KH_{NO2}+S_{NO2}})MonodNH4*MonodPO4*X_{H}`"
   "Process 13", ":math:`ρ_13 = µ_{H}hH_{NO}(\frac{KH4_{O2}}{KH4_{O2}+S_{O2}})(\frac{S_{A}}{KH_{A4}+S_{A}})(\frac{S_{A}}{S_{F}+S_{A}})(\frac{S_{NO}}{KH_{NO}+S_{NO}+S_{NO}^{2}/KI_{NO}})MonodNH4*MonodPO4*X_{H}`"
   "Process 14", ":math:`ρ_14 = µ_{H}hH_{N2O}(\frac{KH5_{O2}}{KH5_{O2}+S_{O2}})(\frac{S_{A}}{KH_{A5}+S_{A}})(\frac{S_{A}}{S_{F}+S_{A}})(\frac{S_{N2O}}{KH_{N2O}+S_{N2O}})MonodNH4*MonodPO4*X_{H}`"
   "Process 15", ":math:`ρ_15 = q_{fe}(\frac{KH_{O2}}{KH_{O2}+S_{O2}})(\frac{KH_{NO3}}{KH_{NO3}+S_{NO3}+S_{NO2}})(\frac{S_{F}}{K_{fe}+S_{F}})X_{H}`"
   "Process 16 (with decay)", ":math:`ρ_16 = b_{H}((/frac{S_{O2}}{KH_{O2}+S_{O2}})+(/frac{hH_{NO3,end}KH_{O2}}{KH_{O2}+(/frac{S_{O2}S_{NO3}{KH_{NO3}+S_{NO3}})})*X_{H}`"
   "Process 16 (without decay)", ":math:`ρ_16 = b_{H}X_{H}`"
   "Process 17", ":math:`ρ_17 = q_{PHA}(/frac{S_{A}}{KP_{A}+S_{A}})MonodPP*X_{PAO}`"
   "Process 18", ":math:`ρ_18 = q_{PP}(/frac{S_{O2}}{KP_{O2}+S_{O2}})MonodSPO4*MonodPHA*MonodIPP*X_{PAO}`"
   "Process 19", ":math:`ρ_19 = q_{PP}*hP_{NO3}(/frac{KP2_{O2}}{KP2_{O2}+S_{O2}})(/frac{S_{NO3}}{KP_{NO3}+S_{NO3}})MonodSPO4_{P}*MonodPHA*MonodIPP*X_{PAO}`"
   "Process 20", ":math:`ρ_20 = q_{PP}*hP_{NO2}(/frac{KP3_{O2}}{KP3_{O2}+S_{O2}})(/frac{S_{NO2}}{KP_{NO2}+S_{NO2}})MonodSPO4_{P}*MonodPHA*MonodIPP*X_{PAO}`"
   "Process 21", ":math:`ρ_21 = q_{PP}*hP_{NO}(/frac{KP4_{O2}}{KP4_{O2}+S_{O2}})(/frac{S_{NO}}{KP_{NO}+S_{NO}+S_{NO}^{2}/KI_{NO}})MonodSPO4_{P}*MonodPHA*MonodIPP*X_{PAO}`"
   "Process 22", ":math:`ρ_22 = q_{PP}*hP_{N2O}(/frac{KP5_{O2}}{KP5_{O2}+S_{O2}})(/frac{S_{N2O}}{KP_{N2O}+S_{N2O}})MonodSPO4_{P}*MonodPHA*MonodIPP*X_{PAO}`"
   "Process 23", ":math:`ρ_23 = µ_{PAO}(/frac{S_{O2}}{KP_{O2}+S_{O2}})MonodSNH4_{P}(/frac{S_{PO4}}{KP_{PO4}+S_{PO4}})MonodPHA*X_{PAO}`"
   "Process 24", ":math:`ρ_24 = µ_{PAO}*hP_{NO3}(/frac{KP2_{O2}}{KP2_{O2}+S_{O2}})(/frac{S_{NO3}}{KP_{NO3}+S_{NO3}})MonodSNH4_{P}*MonodSPO4_{P}*MonodPHA*X_{PAO}`"
   "Process 25", ":math:`ρ_25 = µ_{PAO}*hP_{NO2}(/frac{KP3_{O2}}{KP3_{O2}+S_{O2}})(/frac{S_{NO2}}{KP_{NO2}+S_{NO2}})MonodSNH4_{P}*MonodSPO4_{P}*MonodPHA*X_{PAO}`"
   "Process 26", ":math:`ρ_26 = µ_{PAO}*hP_{NO}(/frac{KP4_{O2}}{KP4_{O2}+S_{O2}})(/frac{S_{NO}}{KP_{NO}+S_{NO}+S_{NO}^{2}/KI_{NO}})MonodSNH4_{P}*MonodSPO4_{P}*MonodPHA*X_{PAO}`"
   "Process 27", ":math:`ρ_27 = µ_{PAO}*hP_{N2O}(/frac{KP5_{O2}}{KP5_{O2}+S_{O2}})(/frac{S_{N2O}}{KP_{N2O}+S_{N2O}})MonodSNH4_{P}*MonodSPO4_{P}*MonodPHA*X_{PAO}`"
   "Process 28 (with decay)", ":math:`ρ_28 = b_{PAO}((/frac{S_{O2}}{KP_{O2}+S_{O2}})+(/frac{hP_{NO3,end}KP_{O2}}{KP_{O2}+(/frac{S_{O2}S_{NO3}{KP_{NO3}+S_{NO3}})})*X_{PAO}`"
   "Process 28 (without decay)", ":math:`ρ_28 = b_{PAO}X_{PAO}`"
   "Process 29 (with decay)", ":math:`ρ_29 = b_{PP}((/frac{S_{O2}}{KP_{O2}+S_{O2}})+(/frac{hPP_{NO3,end}KP_{O2}}{KP_{O2}+(/frac{S_{O2}S_{NO3}{KP_{NO3}+S_{NO3}})})*X_{PP}`"
   "Process 29 (without decay)", ":math:`ρ_29 = b_{PP}X_{Pp}`"
   "Process 30 (with decay)", ":math:`ρ_30 = b_{PHA}((/frac{S_{O2}}{KP_{O2}+S_{O2}})+(/frac{hPHA_{NO3,end}KP_{O2}}{KP_{O2}+(/frac{S_{O2}S_{NO3}{KP_{NO3}+S_{NO3}})})*X_{PHA}`"
   "Process 30 (without decay)", ":math:`ρ_30 = b_{PHA}X_{PHA}`"
   "Process 31", ":math:`ρ_31 = q_{AOB,AMO}(/frac{S_{O2}}{KAOB1_{O2}+S_{O2}})(/frac{S_{FA}}{KAOB_{NH4}+S_{FA}})*X_{AOB}`"
   "Process 32", ":math:`ρ_32 = µ_{AOB,HAO}(/frac{S_{O2}}{KAOB2_{O2}+S_{O2}})(/frac{S_{NH2OH}}{KAOB_{NH2OH}+S_{NH2OH}})(/frac{S_{NH4}}{S_{NH4}+1e-12})(/frac{S_{PO4}}{KAOB_{P}+S_{PO4}})(/frac{X_{AOB}}{1+S_{IS}/K_I_IS_XAOB})`"
   "Process 33", ":math:`ρ_33 = q_{AOB,HAO}(/frac{S_{O2}}{KAOB2_{O2}+S_{O2}})(/frac{S_{NO}}{KAOB_{HAO,NO}+S_{NO}})*X_{AOB}`"
   "Process 34", ":math:`ρ_34 = q_{AOB,N2O,NN}(/frac{S_{NH2OH}}{KAOB_{NH2OH}+S_{NH2OH}})(/frac{S_{NO}}{KAOB_{NN,NO}+S_{NO}})*X_{AOB}`"
   "Process 35", ":math:`ρ_35 = q_{AOB,N2O,ND}(/frac{S_{NH2OH}}{KAOB_{NH2OH}+S_{NH2OH}})(/frac{S_{FNA}}{KAOB_{NH2OH}+S_FNA})*fSO2*X_{AOB}`"
   "Process 36", ":math:`ρ_36 = µ_{NOB}(/frac{S_{O2}}{KNOB_{O2}+S_{O2}})(/frac{S_{FNA}}{KNOB_{NO2}+S_FNA})*(/frac{S_{PO4}}{KNOB_{P}+S_{PO4}})(/frac{X_{NOB}}{1+S_{IS}/K_{I,IS,XNOB}})*X_{AOB}`"
   "Process 37 (with decay)", ":math:`ρ_37 = b_{AOB}((/frac{S_{O2}}{KAOB1_{O2}+S_{O2}})+(/frac{hAOB_{NO3,end}KAOB1_{O2}}{KAOB1_{O2}+S_{O2}})(/frac{S_{NO3}}{KAOB_{NO3}+S_{NO3}})*X_{AOB}`"
   "Process 37 (without decay)", ":math:`ρ_37 = b_{AOB}X_{AOB}`"
   "Process 38 (with decay)", ":math:`ρ_38 = b_{NOB}((/frac{S_{O2}}{KNOB_{O2}+S_{O2}})+(/frac{hNOB_{NO3,end}KNOB_{O2}}{KNOB_{O2}+S_{O2}})(/frac{S_{NO3}}{KNOB_{NO3}+S_{NO3}})*X_{NOB}`"
   "Process 38 (without decay)", ":math:`ρ_38 = b_{NOB}X_{NOB}`"

References
----------
[1] B. Solis, A. Guisasola, X. Flores-Alsina, U. Jeppsson, J.A. Baeza,
A plant-wide model describing GHG emissions and nutrient recovery options for water resource recovery facilities,
Water Research 215 (2022) https://www.sciencedirect.com/science/article/pii/S0043135422001865


