ADM1 to ASM1 Translator
=======================

Introduction
------------

A link is required to translate between biological based and physical or chemical mediated processes
to develop whole-plant modeling of wastewater treatment. This model mediates the interaction between
the Anaerobic Digestor Model 1 (ADM1) and the Activated Sludge Model 1 (ASM1).

The model relies on the following key assumption:

   * supports only liquid phase
   * supports only ADM1 to ASM1 translations

.. index::
   pair: watertap.unit_models.translator_adm1_asm1;translator_adm1_asm1

.. currentmodule:: watertap.unit_models.translator_adm1_asm1

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
.. currentmodule:: watertap.unit_models.translator_adm1_asm1

.. autoclass:: TranslatorData
    :members:
    :noindex:

References
----------
[1] Copp J. and Jeppsson, U., Rosen, C., 2006.
Towards an ASM1 - ADM1 State Variable Interface for Plant-Wide Wastewater Treatment Modeling.
Proceedings of the Water Environment Federation, 2003, pp 498-510.
https://www.accesswater.org/publications/-290550/towards-an-asm1--ndash--adm1-state-variable-interface-for-plant-wide-wastewater-treatment-modeling
