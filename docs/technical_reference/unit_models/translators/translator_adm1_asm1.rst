ADM1 to ASM1 Translator
=======================

Introduction
------------

A link is required to translate between biological and physically- or chemically-mediated processes
to develop whole-plant modeling of wastewater treatment. This model mediates the interaction between
the Anaerobic Digester Model 1 (ADM1) and the Activated Sludge Model 1 (ASM1).

The model relies on the following key assumptions:

   * supports only liquid phase
   * supports only ADM1 to ASM1 translations

.. index::
   pair: watertap.unit_models.translators.translator_adm1_asm1;translator_adm1_asm1

.. currentmodule:: watertap.unit_models.translators.translator_adm1_asm1

Degrees of Freedom
------------------
The translator degrees of freedom are the inlet feed state variables:

    * temperature
    * pressure
    * volumetric flowrate
    * solute compositions
    * cations
    * anions

Ports
-----

This model provides two ports:

* inlet
* outlet

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq']"
   "Inlet Components", ":math:`j`", "['H2O', 'S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_I', 'X_c', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'S_cat', 'S_an', 'S_co2']"
   "Ion", ":math:`j`", "['S_cat', 'S_an'] \  :sup:`*`"
   "Outlet Components", ":math:`j`", "['H2O', 'S_I', 'S_S', 'X_I', 'X_S', 'X_BH', 'X_BA', 'X_P', 'S_O', 'S_NO', 'S_NH', 'S_ND', 'X_ND', 'S_ALK']"
   "Readily Biodegradable COD", ":math:`k`", "['S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac']"
   "Slowly Biodegradable COD", ":math:`m`", "['X_c','X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2']"
        
**Notes**
 :sup:`*` Ion" is a subset of "Component" and uses the same symbol j.

ADM1 Components
---------------
Additional documentation on the ADM1 property model can be found here: `Anaerobic Digestion Model 1 Documentation <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/ADM1.html>`_

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

ASM1 Components
---------------
Additional documentation on the ASM1 property model can be found here: `Activated Sludge Model 1 Documentation <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/ASM1.html>`_

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

.. _Translator_ADM1_ASM1_equations:

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Pressure balance", ":math:`P_{out} = P_{in}`"
   "Temperature balance", ":math:`T_{out} = T_{in}`"
   "Volumetric flow equality", ":math:`F_{out} = F_{in}`"
   "Inert soluble COD balance", ":math:`S_{I, out} = S_{I, in}`"
   "Inert particulate COD balance", ":math:`X_{I, out} = X_{I, in}`"
   "Soluble biodegradable COD balance", ":math:`S_{S, out} = S_{su, in} + S_{aa, in} + S_{fa, in} + S_{va, in} + S_{bu, in} + S_{pro, in} + S_{ac, in}`"
   "Particulate biodegradable COD balance", ":math:`X_{S, out} = X_{c, in} + X_{ch, in} + X_{pr, in} + X_{li, in} + X_{su, in} + X_{aa, in} + X_{fa, in} + X_{c4, in} + X_{pro, in} + X_{ac, in} + X_{h2, in}`"
   "Inorganic nitrogen balance", ":math:`S_{NH, out} = S_{IN, in}`"
   "Organic nitrogen balance", ":math:`S_{ND, out} = (S_{I, in} * N_i) + (S_{aa, in} * N_{aa})`"
   "Particulate organic nitrogen balance", ":math:`X_{ND, out} = N_{bac} * (X_{su, in} + X_{aa, in} + X_{fa, in} + X_{c4, in} + X_{pro, in} + X_{ac, in} + X_{h2, in}) + N_i * X_{I, in} + N_{xc} * X_{c, in} + N_{aa} * X_{pr, in} - x_{ie} * X_{I, in}`"
   "Alkalinity equation", ":math:`S_{ALK, out} = S_{IC, in}`"

Classes
-------
.. currentmodule:: watertap.unit_models.translators.translator_adm1_asm1

.. autoclass:: TranslatorDataADM1ASM1
    :members:
    :noindex:

References
----------
[1] Copp J. and Jeppsson, U., Rosen, C., 2006.
Towards an ASM1 - ADM1 State Variable Interface for Plant-Wide Wastewater Treatment Modeling.
Proceedings of the Water Environment Federation, 2003, pp 498-510.
https://www.accesswater.org/publications/proceedings/-290550/towards-an-asm1---adm1-state-variable-interface-for-plant-wide-wastewater-treatment-modeling
