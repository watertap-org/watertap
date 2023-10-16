ADM1 Property Package
=====================

This package implements properties and reactions of an anaerobic digestion model for industrial wastewater treatment from food and beverage production using an anaerobic digester as provided in
`Batstone, D. J. et al. (2002) <https://iwaponline.com/wst/article-abstract/45/10/65/6034>`_ and `Rosen and Jeppsson (2006) <https://www.iea.lth.se/WWTmodels_download/TR_ADM1.pdf>`_.

This Anaerobic Digestion Model no.1 (ADM1) property/reaction package:
   * supports 'H2O', 'S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_I', 'X_c', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'S_cat', 'S_an', and 'S_co2' as components
   * supports only liquid and vapor phase

Sets
----
.. csv-table::
  :header: "Description", "Symbol", "Indices"

  "Components", ":math:`j`", "['H2O', 'S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_I', 'X_c', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'S_cat', 'S_an', 'S_co2']"
  "Phases", ":math:`p`", "['Liq', 'Vap']"

Components
----------
.. csv-table::
  :header: "Description", "Symbol", "Variable"

  "Monosaccharides, S_su", ":math:`S_{su}`", "S_su"
  "Amino acids, S_aa", ":math:`S_{aa}`", "S_aa"
  "Long chain fatty acids, S_fa", ":math:`S_{fa}`", "S_fa"
  "Total valerate, S_va", ":math:`S_{va}`", "S_va"
  "Total butyrate, S_bu", ":math:`S_{bu}`", "S_bu"
  "Total propionate, S_pro", ":math:`S_{pro}`", "S_pro"
  "Total acetate, S_ac", ":math:`S_{ac}`", "S_ac"
  "Hydrogen gas, S_h2", ":math:`S_{h2}`", "S_h2"
  "Methane gas, S_ch4", ":math:`S_{ch4}`", "S_ch4"
  "Inorganic carbon, S_IC", ":math:`S_{IC}`", "S_IC"
  "Inorganic nitrogen, S_IN", ":math:`S_{IN}`", "S_IN"
  "Soluble inerts, S_I", ":math:`S_I`", "S_I"
  "Composites, X_c", ":math:`X_c`", "X_c"
  "Carbohydrates, X_ch", ":math:`X_{ch}`", "X_ch"
  "Proteins, X_pr", ":math:`X_{pr}`", "X_pr"
  "Lipids, X_li", ":math:`X_{li}`", "X_li"
  "Sugar degraders, X_su", ":math:`X_{su}`", "X_su"
  "Amino acid degraders, X_aa", ":math:`X_{aa}`", "X_aa"
  "Long chain fatty acid (LCFA) degraders, X_fa", ":math:`X_{fa}`", "X_fa"
  "Valerate and butyrate degraders, X_c4", ":math:`X_{c4}`", "X_c4"
  "Propionate degraders, X_pro", ":math:`X_{pro}`", "X_pro"
  "Acetate degraders, X_ac", ":math:`X_{ac}`", "X_ac"
  "Hydrogen degraders, X_h2", ":math:`X_{h2}`", "X_h2"
  "Particulate inerts, X_I", ":math:`X_I`", "X_I"
  "Total cation equivalents concentration, S_cat", ":math:`S_{cat}`", "S_cat"
  "Total anion equivalents concentration, S_an", ":math:`S_{an}`", "S_an"
  "Carbon dioxide, S_co2", ":math:`S_{co2}`", "S_co2"

**NOTE: S_h2 and S_ch4 have vapor phase and liquid phase, S_co2 only has vapor phase, and the other components only have liquid phase. The amount of CO2 dissolved in the liquid phase is equivalent to S_IC - S_HCO3-.**

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Total volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"
   "Component mass concentrations", ":math:`C_j`", "conc_mass_comp", "[j]", ":math:`\text{kg/}\text{m}^3`"
   "Anions in molar concentrations", ":math:`M_a`", "anions", "None", ":math:`\text{kmol/}\text{m}^3`"
   "Cations in molar concentrations", ":math:`M_c`", "cations", "None", ":math:`\text{kmol/}\text{m}^3`"
   "Component pressure", ":math:`P_{j,sat}`", "pressure_sat", "[j]", ":math:`\text{Pa}`"
   "Reference temperature", ":math:`T_{ref}`", "temperature_ref", "None", ":math:`\text{K}`"
   "Reference component mass concentrations", ":math:`C_{j,ref}`", "conc_mass_comp_ref", "[j]", ":math:`\text{kg/}\text{m}^3`"

Stoichiometric Parameters
-------------------------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Value at 20 C", "Units"

   "Soluble inerts from composites, f_sI_xc", ":math:`f_{sI,xc}`", "f_sI_xc", 0.1, ":math:`\text{dimensionless}`"
   "Particulate inerts from composites, f_xI_xc", ":math:`f_{xI,xc}`", "f_xI_xc", 0.2, ":math:`\text{dimensionless}`"
   "Carbohydrates from composites, f_ch_xc", ":math:`f_{ch,xc}`", "f_ch_xc", 0.2, ":math:`\text{dimensionless}`"
   "Proteins from composites, f_pr_xc", ":math:`f_{pr,xc}`", "f_pr_xc", 0.2, ":math:`\text{dimensionless}`"
   "Lipids from composites, f_li_xc", ":math:`f_{li,xc}`", "f_li_xc", 0.3, ":math:`\text{dimensionless}`"
   "Nitrogen content of composites, N_xc", ":math:`N_{xc}`", "N_xc", 0.0376/14, ":math:`\text{kmol-N/}\text{kg-COD}`"
   "Nitrogen content of inerts, N_I", ":math:`N_I`", "N_I", 0.06/14, ":math:`\text{kmol-N/}\text{kg-COD}`"
   "Nitrogen in amino acids and proteins, N_aa", ":math:`N_{aa}`", "N_aa", 0.007, ":math:`\text{kmol-N/}\text{kg-COD}`"
   "Nitrogen content in bacteria, N_bac", ":math:`N_{bac}`", "N_bac", 0.08/14, ":math:`\text{kmol-N/}\text{kg-COD}`"
   "Fatty acids from lipids, f_fa_li", ":math:`f_{fa,li}`", "f_fa_li", 0.95, ":math:`\text{dimensionless}`"
   "Hydrogen from sugars, f_h2_su", ":math:`f_{h2,su}`", "f_h2_su", 0.19, ":math:`\text{dimensionless}`"
   "Butyrate from sugars, f_bu_su", ":math:`f_{bu,su}`", "f_bu_su", 0.13, ":math:`\text{dimensionless}`"
   "Propionate from sugars, f_pro_su", ":math:`f_{pro,su}`", "f_pro_su", 0.27, ":math:`\text{dimensionless}`"
   "Acetate from sugars, f_ac_su", ":math:`f_{ac,su}`", "f_ac_su", 0.41, ":math:`\text{dimensionless}`"
   "Hydrogen from amino acids, f_h2_aa", ":math:`f_{h2,aa}`", "f_h2_aa", 0.06, ":math:`\text{dimensionless}`"
   "Valerate from amino acids, f_va_aa", ":math:`f_{va,aa}`", "f_va_aa", 0.23, ":math:`\text{dimensionless}`"
   "Butyrate from amino acids, f_bu_aa", ":math:`f_{bu,aa}`", "f_bu_aa", 0.26, ":math:`\text{dimensionless}`"
   "Propionate from amino acids, f_pro_aa", ":math:`f_{pro,aa}`", "f_pro_aa", 0.05, ":math:`\text{dimensionless}`"
   "Acetate from amino acids, f_ac_aa", ":math:`f_{ac,aa}`", "f_ac_aa", 0.4, ":math:`\text{dimensionless}`"
   "Yield of biomass on sugar substrate, Y_su", ":math:`Y_{su}`", "Y_su", 0.1, ":math:`\text{kg-COD X/}\text{kg-COD S}`"
   "Yield of biomass on amino acid substrate, Y_aa", ":math:`Y_{aa}`", "Y_aa", 0.08, ":math:`\text{kg-COD X/}\text{kg-COD S}`"
   "Yield of biomass on fatty acid substrate, Y_fa", ":math:`Y_{fa}`", "Y_fa", 0.06, ":math:`\text{kg-COD X/}\text{kg-COD S}`"
   "Yield of biomass on valerate and butyrate substrate, Y_c4", ":math:`Y_{c4}`", "Y_c4", 0.06, ":math:`\text{kg-COD X/}\text{kg-COD S}`"
   "Yield of biomass on propionate substrate, Y_pro", ":math:`Y_{pro}`", "Y_pro", 0.04, ":math:`\text{kg-COD X/}\text{kg-COD S}`"
   "Yield of biomass on acetate substrate, Y_ac", ":math:`Y_{ac}`", "Y_ac", 0.05, ":math:`\text{kg-COD X/}\text{kg-COD S}`"
   "Yield of hydrogen per biomass, Y_h2", ":math:`Y_{h2}`", "Y_h2", 0.06, ":math:`\text{kg-COD X/}\text{kg-COD S}`"

Kinetic Parameters
------------------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Value at 20 C", "Units"

   "First-order kinetic parameter for disintegration, k_dis", ":math:`k_{dis}`", "k_dis", 0.5, ":math:`\text{d}^{-1}`"
   "First-order kinetic parameter for hydrolysis of carbohydrates, k_hyd_ch", ":math:`k_{hyd,ch}`", "k_hyd_ch", 10, ":math:`\text{d}^{-1}`"
   "First-order kinetic parameter for hydrolysis of proteins, k_hyd_pr", ":math:`k_{hyd,pr}`", "k_hyd_pr", 10, ":math:`\text{d}^{-1}`"
   "First-order kinetic parameter for hydrolysis of lipids, k_hyd_li", ":math:`k_{hyd,li}`", "k_hyd_li", 10, ":math:`\text{d}^{-1}`"
   "Inhibition parameter for inorganic nitrogen, K_S_IN", ":math:`K_{S_{IN}}`", "K_S_IN", 1e-4, ":math:`\text{kmol/}\text{m}^3`"
   "Monod maximum specific uptake rate of sugars, k_m_su", ":math:`k_{m_{su}}`", "k_m_su", 30, ":math:`\text{d}^{-1}`"
   "Half saturation value for uptake of sugars, K_S_su", ":math:`K_{S_{su}}`", "K_S_su", 0.5, ":math:`\text{kg/}\text{m}^3`"
   "Upper limit of pH for uptake rate of amino acids, pH_UL_aa", ":math:`pH_{UL,aa}`", "pH_UL_aa", 5.5, ":math:`\text{dimensionless}`"
   "Lower limit of pH for uptake rate of amino acids, pH_LL_aa", ":math:`pH_{LL,aa}`", "pH_LL_aa", 4, ":math:`\text{dimensionless}`"
   "Monod maximum specific uptake rate of amino acids, k_m_aa", ":math:`k_{m_{aa}}`", "k_m_aa", 50, ":math:`\text{d}^{-1}`"
   "Half saturation value for uptake of amino acids, K_S_aa", ":math:`K_{S_{aa}}`", "K_S_aa", 0.3, ":math:`\text{kg/}\text{m}^3`"
   "Monod maximum specific uptake rate of fatty acids, k_m_fa", ":math:`k_{m_{fa}}`", "k_m_fa", 6, ":math:`\text{d}^{-1}`"
   "Half saturation value for uptake of fatty acids, K_S_fa", ":math:`K_{S_{fa}}`", "K_S_fa", 0.4, ":math:`\text{kg/}\text{m}^3`"
   "Inhibition parameter for hydrogen during uptake of fatty acids, K_I_h2_fa", ":math:`K_{I,h2_{fa}}`", "K_I_h2_fa", 5e-6, ":math:`\text{kg/}\text{m}^3`"
   "Monod maximum specific uptake rate of valerate and butyrate, k_m_c4", ":math:`k_{m_{c4}}`", "k_m_c4", 20, ":math:`\text{d}^{-1}`"
   "Half saturation value for uptake of valerate and butyrate, K_S_c4", ":math:`K_{S_{c4}}`", "K_S_c4", 0.2, ":math:`\text{kg/}\text{m}^3`"
   "Inhibition parameter for hydrogen during uptake of valerate and butyrate, K_I_h2_c4", ":math:`K_{I,h2_{c4}}`", "K_I_h2_c4", 1e-5, ":math:`\text{kg/}\text{m}^3`"
   "Monod maximum specific uptake rate of propionate, k_m_pro", ":math:`k_{m_{pro}}`", "k_m_pro", 13, ":math:`\text{d}^{-1}`"
   "Half saturation value for uptake of propionate, K_S_pro", ":math:`K_{S_{pro}}`", "K_S_pro", 0.1, ":math:`\text{kg/}\text{m}^3`"
   "Inhibition parameter for hydrogen during uptake of propionate, K_I_h2_pro", ":math:`K_{I,h2_{pro}}`", "K_I_h2_pro", 3.5e-6, ":math:`\text{kg/}\text{m}^3`"
   "Monod maximum specific uptake rate of acetate, k_m_ac", ":math:`k_{m_{ac}}`", "k_m_ac", 8, ":math:`\text{d}^{-1}`"
   "Half saturation value for uptake of acetate, K_S_ac", ":math:`K_{S_{ac}}`", "K_S_ac", 0.15, ":math:`\text{kg/}\text{m}^3`"
   "Inhibition parameter for ammonia during uptake of acetate, K_I_nh3", ":math:`K_{I,nh3}`", "K_I_nh3", 0.0018, ":math:`\text{kg/}\text{m}^3`"
   "Upper limit of pH for uptake rate of acetate, pH_UL_ac", ":math:`pH_{UL,ac}`", "pH_UL_ac", 7, ":math:`\text{dimensionless}`"
   "Lower limit of pH for uptake rate of acetate, pH_LL_ac", ":math:`pH_{LL,ac}`", "pH_LL_ac", 6, ":math:`\text{dimensionless}`"
   "Monod maximum specific uptake rate of hydrogen, k_m_h2", ":math:`k_{m_{h2}}`", "k_m_h2", 35, ":math:`\text{d}^{-1}`"
   "Half saturation value for uptake of hydrogen, K_S_h2", ":math:`K_{S_{h2}}`", "K_S_h2", 7e-6, ":math:`\text{kg/}\text{m}^3`"
   "Upper limit of pH for uptake rate of hydrogen, pH_UL_h2", ":math:`pH_{UL,h2}`", "pH_UL_h2", 6, ":math:`\text{dimensionless}`"
   "Lower limit of pH for uptake rate of hydrogen, pH_LL_h2", ":math:`pH_{LL,h2}`", "pH_LL_h2", 5, ":math:`\text{dimensionless}`"
   "First-order decay rate for X_su, k_dec_X_su", ":math:`k_{dec,X_{su}}`", "k_dec_X_su", 0.02, ":math:`\text{d}^{-1}`"
   "First-order decay rate for X_aa, k_dec_X_aa", ":math:`k_{dec,X_{aa}}`", "k_dec_X_aa", 0.02, ":math:`\text{d}^{-1}`"
   "First-order decay rate for X_fa, k_dec_X_fa", ":math:`k_{dec,X_{fa}}`", "k_dec_X_fa", 0.02, ":math:`\text{d}^{-1}`"
   "First-order decay rate for X_c4, k_dec_X_c4", ":math:`k_{dec,X_{c4}}`", "k_dec_X_c4", 0.02, ":math:`\text{d}^{-1}`"
   "First-order decay rate for X_pro, k_dec_X_pro", ":math:`k_{dec,X_{pro}}`", "k_dec_X_pro", 0.02, ":math:`\text{d}^{-1}`"
   "First-order decay rate for X_ac, k_dec_X_ac", ":math:`k_{dec,X_{ac}}`", "k_dec_X_ac", 0.02, ":math:`\text{d}^{-1}`"
   "First-order decay rate for X_h2, k_dec_X_h2", ":math:`k_{dec,X_{h2}}`", "k_dec_X_h2", 0.02, ":math:`\text{d}^{-1}`"
   "Water dissociation constant, pK_W", ":math:`pK_W`", "pKW", 14, ":math:`\text{dimensionless}`"
   "Process inhibition term, I", ":math:`I`", "I", 1, ":math:`\text{dimensionless}`"
   "Valerate acid-base equilibrium constant, K_a_va", ":math:`K_{a,va}`", "K_a_va", 1.38e-5, ":math:`\text{kmol/}\text{m}^3`"
   "Butyrate acid-base equilibrium constant, K_a_bu", ":math:`K_{a,bu}`", "K_a_bu", 1.5e-5, ":math:`\text{kmol/}\text{m}^3`"
   "Propionate acid-base equilibrium constant, K_a_pro", ":math:`K_{a,pro}`", "K_a_bu", 1.32e-5, ":math:`\text{kmol/}\text{m}^3`"
   "Acetate acid-base equilibrium constant, K_a_ac", ":math:`K_{a,ac}`", "K_a_ac", 1.74e-5, ":math:`\text{kmol/}\text{m}^3`"
   "Carbon dioxide acid-base equilibrium constant, pK_a_co2", ":math:`pK_{a,co2}`", "pK_a_co2", 6.35, ":math:`\text{dimensionless}`"
   "Inorganic nitrogen acid-base equilibrium constant, pK_a_IN", ":math:`pK_{a,IN}`", "pK_a_IN", 9.25, ":math:`\text{dimensionless}`"

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

   "Disintegration", ":math:`\rho_1 = k_{dis} C_{X_c}`"
   "Hydrolysis of carbohydrates", ":math:`\rho_2 = k_{hyd,ch} C_{X_{ch}}`"
   "Hydrolysis of proteins", ":math:`\rho_3 = k_{hyd,pr} C_{X_{pr}}`"
   "Hydrolysis of lipids", ":math:`\rho_4 = k_{hyd,li} C_{X_{li}}`"
   "Uptake of sugars", ":math:`\rho_5 = k_{m_{su}} \frac{C_{S_{su}}}{K_{S_{su}}+C_{S_{su}}} C_{X_{su}} \frac{1}{1 + K_{S_{IN}}/C_{S_{IN}}/14} I_{pH,aa}`"
   "Uptake of amino acids", ":math:`\rho_6 = k_{m_{aa}} \frac{C_{S_{aa}}}{K_{S_{aa}}+C_{S_{aa}}} C_{X_{aa}} \frac{1}{1 + K_{S_{IN}}/C_{S_{IN}}/14} I_{pH,aa}`"
   "Uptake of long chain fatty acids (LCFAs)", ":math:`\rho_7 = k_{m_{fa}} \frac{C_{S_{fa}}}{K_{S_{fa}}+C_{S_{fa}}} C_{X_{fa}} \frac{1}{1 + K_{S_{IN}}/C_{S_{IN}}/14} \cdot \frac{1}{1 + C_{S_{h2}}/K_{I,h2_{fa}}} I_{pH,aa}`"
   "Uptake of valerate", ":math:`\rho_8 = k_{m_{c4}} \frac{C_{S_{va}}}{K_{S_{c4}}+C_{S_{va}}} C_{X_{c4}} \frac{C_{S_{va}}}{C_{S_{bu}} + C_{S_{va}}} \cdot \frac{1}{1 + K_{S_{IN}}/C_{S_{IN}}/14} \cdot \frac{1}{1 + C_{S_{h2}}/K_{I,h2_{c4}}} I_{pH,aa}`"
   "Uptake of butyrate", ":math:`\rho_9 = k_{m_{c4}} \frac{C_{S_{bu}}}{K_{S_{c4}}+C_{S_{bu}}} C_{X_{c4}} \frac{C_{S_{bu}}}{C_{S_{bu}} + C_{S_{va}}} \cdot \frac{1}{1 + K_{S_{IN}}/C_{S_{IN}}/14} \cdot \frac{1}{1 + C_{S_{h2}}/K_{I,h2_{c4}}} I_{pH,aa}`"
   "Uptake of propionate", ":math:`\rho_{10} = k_{m_{pro}} \frac{C_{S_{pro}}}{K_{S_{pro}}+C_{S_{pro}}} C_{X_{pro}} \cdot \frac{1}{1 + K_{S_{IN}}/C_{S_{IN}}/14} \cdot \frac{1}{1 + C_{S_{h2}}/K_{I,h2_{pro}}} I_{pH,aa}`"
   "Uptake of acetate", ":math:`\rho_{11} = k_{m_{ac}} \frac{C_{S_{ac}}}{K_{S_{ac}}+C_{S_{ac}}} C_{X_{ac}} \cdot \frac{1}{1 + K_{S_{IN}}/C_{S_{IN}}/14} \cdot \frac{1}{1 + C_{NH3}/K_{I,nh3}} I_{pH,ac}`"
   "Uptake of hydrogen", ":math:`\rho_{12} = k_{m_{h2}} \frac{C_{S_{h2}}}{K_{S_{h2}}+C_{S_{h2}}} C_{X_{h2}} \cdot \frac{1}{1 + K_{S_{IN}}/C_{S_{IN}}/14} I_{pH,h2}`"
   "Decay of X_su", ":math:`\rho_{13} = k_{dec, X_{su}} C_{X_{su}}`"
   "Decay of X_aa", ":math:`\rho_{14} = k_{dec, X_{aa}} C_{X_{aa}}`"
   "Decay of X_fa", ":math:`\rho_{15} = k_{dec, X_{fa}} C_{X_{fa}}`"
   "Decay of X_c4", ":math:`\rho_{16} = k_{dec, X_{c4}} C_{X_{c4}}`"
   "Decay of X_pro", ":math:`\rho_{17} = k_{dec, X_{pro}} C_{X_{pro}}`"
   "Decay of X_ac", ":math:`\rho_{18} = k_{dec, X_{ac}} C_{X_{ac}}`"
   "Decay of X_h2", ":math:`\rho_{19} = k_{dec, X_{h2}} C_{X_{h2}}`"

Additional Variables
--------------------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Value at 20 C", "Units"

   "Molar concentration of hydrogen, S_H", ":math:`S_{H}`", "S_H", 1e-7, ":math:`\text{kmol/}\text{m}^3`"
   "pH of solution", ":math:`pH`", "pH", 7, ":math:`\text{dimensionless}`"
   "Mass concentration of valerate, va-", ":math:`C_{va}`", "conc_mass_va", 0.011, ":math:`\text{kg/}\text{m}^3`"
   "Mass concentration of butyrate, bu-", ":math:`C_{bu}`", "conc_mass_bu", 0.013, ":math:`\text{kg/}\text{m}^3`"
   "Mass concentration of propionate, pro-", ":math:`C_{pro}`", "conc_mass_pro", 0.016, ":math:`\text{kg/}\text{m}^3`"
   "Mass concentration of acetate, ac-", ":math:`C_{ac}`", "conc_mass_ac", 0.2, ":math:`\text{kg/}\text{m}^3`"
   "Molar concentration of bicarbonate, HCO3", ":math:`M_{hco3}`", "conc_mol_hco3", 0.14, ":math:`\text{kmol/}\text{m}^3`"
   "Molar concentration of ammonia, NH3", ":math:`M_{nh3}`", "conc_mol_nh3", 0.0041, ":math:`\text{kmol/}\text{m}^3`"
   "Molar concentration of carbon dioxide, CO2", ":math:`M_{co2}`", "conc_mol_co2", 0.0099, ":math:`\text{kmol/}\text{m}^3`"
   "Molar concentration of ammonium, NH4", ":math:`M_{nh4}`", "conc_mol_nh4", 0.1261, ":math:`\text{kmol/}\text{m}^3`"

Additional Constraints
----------------------
.. csv-table::
   :header: "Description", "Equation"

   "Water dissociation constant constraint", ":math:`log(10^{-pK_{W}}) = log(10^{-14}) + (\frac{55900}{R} * (\frac{1}{T_{ref}} - \frac{1}{T}))`"
   "CO2 acid-base equilibrium constraint", ":math:`log(10^{-pK_{a,co2}}) = log(10^{-6.35}) + (\frac{7646}{R} * (\frac{1}{T_{ref}} - \frac{1}{T}))`"
   "Nitrogen acid-base equilibrium constraint", ":math:`log(10^{-pK_{a,IN}}) = log(10^{-9.25}) + (\frac{51965}{R} * (\frac{1}{T_{ref}} - \frac{1}{T}))`"
   "pH of solution", ":math:`pH = -log(S_{H})`"
   "Mass concentration of valerate, va-", ":math:`C_{va,ref} = C_{va} (1 + \frac{S_{H}}{K_{a,va}})`"
   "Mass concentration of butyrate, bu-", ":math:`C_{bu,ref} = C_{bu} (1 + \frac{S_{H}}{K_{a,bu}})`"
   "Mass concentration of propionate, pro-", ":math:`C_{pro,ref} = C_{pro} (1 + \frac{S_{H}}{K_{a,pro}})`"
   "Mass concentration of acetate, ac-", ":math:`C_{ac,ref} = C_{ac} (1 + \frac{S_{H}}{K_{a,ac}})`"
   "Molar concentration of bicarbonate, HCO3", ":math:`pK_{a,co2} = log(M_{co2}) - log(M_{hco3}) + pH`"
   "Molar concentration of ammonia, NH3", ":math:`pK_{a,IN} = log(M_{nh4}) - log(M_{nh3}) + pH`"
   "Molar concentration of carbon dioxide, CO2", ":math:`M_{co2} = \frac{C_{S_{IC},ref}}{12} - M_{hco3}`"
   "Molar concentration of ammonium, NH4+", ":math:`M_{nh4} = \frac{C_{S_{IN},ref}}{14} - M_{nh3}`"
   "Molar concentration of hydrogen, H+", ":math:`S_{H} = M_{hco3} + \frac{C_{ac}}{64} + \frac{C_{pro}}{112} + \frac{C_{bu}}{160} + \frac{C_{va}}{208} + 10^{pH - pK_{W}} + M_{a} - M_{c} - M_{nh4}`"

The rules for pH inhibition of amino-acid-utilizing microorganisms (:math:`I_{pH,aa}`), acetate-utilizing microorganisms (:math:`I_{pH,ac}`), and hydrogen-utilizing microorganisms (:math:`I_{pH,h2}`) are:

    .. math::

       I_{pH,aa}=
       \begin{cases}
         \exp{(-3 (\frac{pH - pH_{UL,aa}}{pH_{UL,aa} - pH_{LL,aa}})^2)} & \text{for } pH \le pH_{UL,aa}\\
         1 & \text{for } pH > pH_{UL,aa}
       \end{cases}

       I_{pH,ac}=
       \begin{cases}
         \exp{(-3 (\frac{pH - pH_{UL,ac}}{pH_{UL,ac} - pH_{LL,ac}})^2)} & \text{for } pH \le pH_{UL,ac}\\
         1 & \text{for } pH > pH_{UL,ac}
       \end{cases}

       I_{pH,h2}=
       \begin{cases}
         \exp{(-3 (\frac{pH - pH_{UL,h2}}{pH_{UL,h2} - pH_{LL,h2}})^2)} & \text{for } pH \le pH_{UL,h2}\\
         1 & \text{for } pH > pH_{UL,h2}
       \end{cases}

The rules for inhibition related to secondary substrate (:math:`I_{IN,lim}`), hydrogen inhibition attributed to long chain fatty acids (:math:`I_{h2,fa}`), hydrogen inhibition attributed to valerate and butyrate uptake (:math:`I_{h2,c4}`), hydrogen inhibition attributed to propionate uptake (:math:`I_{h2,pro}`), ammonia inibition attributed to acetate uptake (:math:`I_{nh3}`), are:

    .. math::

       I_{IN,lim} = \frac{1}{1 + \frac{K_{S_{IN}}}{C_{S_{IN}}/14}}

       I_{h2, fa}= \frac{1}{1 + \frac{C_{S_{h2}}}{K_{I,h2,fa}}}

       I_{h2, c4}= \frac{1}{1 + \frac{C_{S_{h2}}}{K_{I,h2,c4}}}

       I_{h2, pro}= \frac{1}{1 + \frac{C_{S_{h2}}}{K_{I,h2,pro}}}

       I_{nh3}= \frac{1}{1 + \frac{M_{nh3}}{K_{I,nh3}}}


Class Documentation
-------------------
.. currentmodule:: watertap.property_models.anaerobic_digestion.adm1_properties

.. autoclass:: ADM1ParameterBlock
    :members:
    :noindex:

.. autoclass:: ADM1ParameterData
    :members:
    :noindex:

.. autoclass:: _ADM1StateBlock
    :members:
    :noindex:

.. autoclass:: ADM1StateBlockData
    :members:
    :noindex:

.. currentmodule:: watertap.property_models.anaerobic_digestion.adm1_properties_vapor

.. autoclass:: ADM1_vaporParameterBlock
    :members:
    :noindex:

.. autoclass:: ADM1_vaporParameterData
    :members:
    :noindex:

.. autoclass:: _ADM1_vaporStateBlock
    :members:
    :noindex:

.. autoclass:: ADM1_vaporStateBlockData
    :members:
    :noindex:

.. currentmodule:: watertap.property_models.anaerobic_digestion.adm1_reactions

.. autoclass:: ADM1ReactionParameterBlock
    :members:
    :noindex:

.. autoclass:: ADM1ReactionParameterData
    :members:
    :noindex:

.. autoclass:: _ADM1ReactionBlock
    :members:
    :noindex:

.. autoclass:: ADM1ReactionBlockData
    :members:
    :noindex:


References
----------
[1] Batstone, D.J., Keller, J., Angelidaki, I., Kalyuzhnyi, S.V., Pavlostathis, S.G., Rozzi, A., Sanders, W.T.M., Siegrist, H.A. and Vavilin, V.A., 2002.
The IWA anaerobic digestion model no 1 (ADM1).
Water Science and technology, 45(10), pp.65-73.
https://iwaponline.com/wst/article-abstract/45/10/65/6034

[2] Rosen, C. and Jeppsson, U., 2006.
Aspects on ADM1 Implementation within the BSM2 Framework.
Department of Industrial Electrical Engineering and Automation, Lund University, Lund, Sweden, pp.1-35.
https://www.iea.lth.se/WWTmodels_download/TR_ADM1.pdf