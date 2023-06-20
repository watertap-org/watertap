ADM1 to ASM2d Translator
========================

Introduction
------------

A link is required to translate between biological based and physical or chemical mediated processes
to develop plant-wide modeling of wastewater treatment. This model mediates the interaction between
the Modified Anaerobic Digestor Model 1 (ADM1) and the Modified Activated Sludge Model 2d (ASM2d).

The model relies on the following key assumption:

   * supports only liquid phase
   * supports only Modified ADM1 to Modified ASM2d translations

.. index::
   pair: watertap.unit_models.translators.translator_adm1_asm2d;translator_adm1_asm2d

.. currentmodule:: watertap.unit_models.translators.translator_adm1_asm2d

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
   "Inlet Components", ":math:`j`", "['H2O', 'S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_IP', 'S_I', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'X_PHA', 'X_PP', 'X_PAO', 'S_K', 'S_Mg']"
   "Ion", ":math:`j`", "['S_cat', 'S_an'] \  :sup:`1`"
   "Outlet Components", ":math:`j`", "['H2O', 'S_A', 'S_F', 'S_I', 'S_N2', 'S_NH4', 'S_NO3', 'S_O2', 'S_PO4', 'S_K', 'S_Mg', 'S_IC', 'X_AUT', 'X_H', 'X_I', 'X_PAO', 'X_PHA', 'X_PP', 'X_S']
   "Readily Biodegradable COD", ":math:`k`", "['S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac']"
   "Slowly Biodegradable COD", ":math:`m`", "['X_ch', 'X_pr', 'X_li']"
   "Unchanged Components", ":math:`j`", "['S_I', 'X_I', 'X_PP', 'X_PHA', 'S_K', 'S_Mg', 'S_IC']"
   "Zero Flow Components", ":math:`j`", "['S_N2', 'S_NO3', 'S_O2', 'X_AUT', 'X_H', 'X_PAO']"

**Notes**
 :sup:`1` Ion" is a subset of "Inlet Components" and uses the same symbol j.

.. _Translator_ADM1_ASM2d_equations:

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Volumetric flow equality", ":math:`F_{out} = F_{in}`"
   "Temperature balance", ":math:`T_{out} = T_{in}`"
   "Pressure balance", ":math:`P_{out} = P_{in}`"
   "Fermentable substrate conversion", ":math:`S_{F, out} = S_{su, in} + S_{aa, in} + S_{fa, in}`"
   "Acetic acid conversion", ":math:`S_{A, out} = S_{va, in} + S_{bu, in} + S_{pro, in} + S_{ac, in}`"
   "Inert soluble COD conversion", ":math:`S_{I, out} = S_{I, in}`"
   "Inert particulate COD conversion", ":math:`X_{I, out} = X_{I, in}`"
   "Polyphosphate conversion", ":math:`X_{PP, out} = X_{PP, in}`"
   "Polyhydroxyalkanoates conversion", ":math:`X_{PHA, out} = X_{PHA, in}`"
   "Potassium conversion", ":math:`S_{K, out} = S_{K, in}`"
   "Magnesium conversion", ":math:`S_{Mg, out} = S_{Mg, in}`"
   "Inorganic carbon conversion", ":math:`S_{IC, out} = S_{IC, in}`"
   "Ammonium conversion", ":math:`S_{NH4, out} = S_{IN, in}`"
   "Phosphate conversion", ":math:`S_{PO4, out} = S_{IP, in}`"
   "Biodegradable particulate organics conversion", ":math:`X_{S, out} = X_{ch, in} + X_{pr, in} + X_{li, in}`"


Classes
-------
.. currentmodule:: watertap.unit_models.translators.translator_adm1_asm2d

.. autoclass:: TranslatorDataADM1ASM2D
    :members:
    :noindex:

References
----------
[1] Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J., 2016.
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes.
Water Research, 95, pp.370-382.

