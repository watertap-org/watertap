ADM1 to ASM2d Translator
========================

.. note::

    Documentation for this model is undergoing refinement.

Introduction
------------

A link is required to translate between biological based and physical or chemical mediated processes
to develop plant-wide modeling of wastewater treatment. This model mediates the interaction between
the Anaerobic Digestor Model 1 (ADM1) and the Activated Sludge Model 2d (ASM2d).

The model relies on the following key assumption:

   * supports only liquid phase
   * the inlet property package is adm1_properties
   * the outlet property package is asm2d_properties

.. index::
   pair: watertap.unit_models.translator_adm1_asm2d;translator_adm1_asm2d

.. currentmodule:: watertap.unit_models.translator_adm1_asm2d

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
   "Inlet Components", ":math:`j`", "['H2O', 'S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac','S_I','S_IN','S_IP','S_IC','X_I','X_ch','X_pr','X_li','X_PP','X_PHA',]"
   "Ion", ":math:`j`", "['S_cat', 'S_an'] \  :sup:`*`"
   "Outlet Components", ":math:`j`", "['H2O', 'S_I','S_F','S_A','S_I','S_NH4','S_PO4','S_IC','X_I','X_S','X_PP','X_PHA',]"
   "Readily Biodegradable COD", ":math:`k`", "['S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac']"
   "Slowly Biodegradable COD", ":math:`m`", "['X_ch', 'X_pr', 'X_li']"
   "Zero Flow Components", ":math:`j`", "['S_N2','S_NO3','S_O2','S_PO4','S_NO','S_ALK','X_AUT','X_H','X_MeOH','X_MeP','X_PAO','X_PHA','X_PP','X_TSS']

**Notes**
 :sup:`*` Ion" is a subset of "Inlet Components" and uses the same symbol j.

.. _Translator_ADM1_ASM2d_equations:

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Pressure balance", ":math:`P_{out} = P_{in}`"
   "Temperature balance", ":math:`T_{out} = T_{in}`"
   "Volumetric flow equality", ":math:`F_{out} = F_{in}`"
   "Fermentable substrate conversion", ":math:`S_{F, out} = S_{su, in} + S_{aa, in} + S_{fa, in}`"
   "Acetic acid conversion", ":math:`S_{A, out} = S_{va, in} + S_{bu, in} + S_{pro, in} + S_{ac, in}`"
   "Inert soluble COD conversion", ":math:`S_{I, out} = S_{I, in}`"
   "Ammonium conversion", ":math:`S_{NH4, out} = S_{IN, in}`"
   "Phosphate conversion", ":math:`S_{PO4, out} = S_{IP, in}`"
   "Inorganic carbon conversion", ":math:`S_{IC, out} = S_{IC, in}`"
   "Inert particulate COD conversion", ":math:`X_{I, out} = X_{I, in}`"
   "Biodegradable particulate organics conversion", ":math:`X_{S, out} = X_{ch, in} + X_{pr, in} + X_{li, in}`"



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
.. currentmodule:: watertap.unit_models.translator_adm1_asm2d

.. autoclass:: UnitModelBlockData
    :members:
    :noindex:

References
----------
[1] Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J., 2016.
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes.
Water Research, 95, pp.370-382.

