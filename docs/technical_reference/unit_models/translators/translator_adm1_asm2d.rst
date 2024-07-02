ADM1 to ASM2d Translator
========================

Introduction
------------

A link is required to translate between biological based and physical or chemical mediated processes
to develop plant-wide modeling of wastewater treatment. This model mediates the interaction between
the Modified Anaerobic Digester Model 1 (ADM1) and the Modified Activated Sludge Model 2d (ASM2d).

The model relies on the following key assumptions:

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
   "Inlet Components", ":math:`j_{in}`", "['H2O', 'S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_IP', 'S_I', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'X_PHA', 'X_PP', 'X_PAO', 'S_K', 'S_Mg']"
   "Outlet Components", ":math:`j_{out}`", "['H2O', 'S_A', 'S_F', 'S_I', 'S_N2', 'S_NH4', 'S_NO3', 'S_O2', 'S_PO4', 'S_K', 'S_Mg', 'S_IC', 'X_AUT', 'X_H', 'X_I', 'X_PAO', 'X_PHA', 'X_PP', 'X_S']"
   "Ion", ":math:`j_{in}`", "['S_cat', 'S_an'] \  :sup:`1`"
   "Zero Flow Components", ":math:`z`", "['S_O2', 'S_N2', 'S_NO3', 'X_H', 'X_PAO', 'X_PP', 'X_PHA', 'X_AUT']"

**Notes**
 :sup:`1` "Ion" is a subset of "Inlet Components" and uses the same symbol j_in.

Parameters
----------

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Value", "Units"

   "Soluble inerts from composites", ":math:`f_{sI, xc}`", "f_sI_xc", 0, ":math:`\text{dimensionless}`"
   "Particulate inerts from composites", ":math:`f_{xI, xc}`", "f_xI_xc", 0.1, ":math:`\text{dimensionless}`"
   "Carbohydrates from composites", ":math:`f_{ch, xc}`", "f_ch_xc", 0.275, ":math:`\text{dimensionless}`"
   "Proteins from composites", ":math:`f_{pr, xc}`", "f_pr_xc", 0.275, ":math:`\text{dimensionless}`"
   "Lipids from composites", ":math:`f_{li, xc}`", "f_li_xc", 0.35, ":math:`\text{dimensionless}`"
   "Phosphorus content of X_ch", ":math:`P_{ch}`", "P_ch", 0, ":math:`\text{dimensionless}`"

.. _Translator_ADM1_ASM2d_equations:

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Volumetric flow equality", ":math:`F_{out} = F_{in}`"
   "Temperature balance", ":math:`T_{out} = T_{in}`"
   "Pressure balance", ":math:`P_{out} = P_{in}`"
   "Zero-flow component conversions", ":math:`C_{z, out} = 0`"
   "Biomass concentration", ":math:`X_{bio} = X_{su, in} + X_{aa, in} + X_{fa, in} + X_{c4, in} + X_{pro, in} + X_{ac, in} + X_{h2, in} + X_{PAO, in}`"
   "S_ac concentration", ":math:`S_{ac, 1} = S_{ac, in} + X_{PHA, in}`"
   "S_IC concentration", ":math:`S_{IC, 1} = (X_{bio} * Ci[X_{su}] * 12) - (X_{bio} * f_{sI, xc} * Ci[S_{I}] * 12) - (X_{bio} * f_{ch, xc} * Ci[X_{ch}] * 12) - (X_{bio} * f_{pr, xc} * Ci[X_{pr}] * 12) - (X_{bio} * f_{li, xc} * Ci[X_{li}] * 12) - (X_{bio} * f_{xI, xc} * Ci[X_{I}] * 12) + (X_{PHA, in} * Ci[X_{PHA}] * 12) - (X_{PHA, in} * Ci[S_{ac}] * 12)`"
   "S_IN concentration", ":math:`S_{IN, 1} = (X_{bio} * Ci[X_{su}] * 14) - (X_{bio} * f_{sI, xc} * Ni[S_{I}] * 14) - (X_{bio} * f_{pr, xc} * Ni[X_{pr}] * 14) - (X_{bio} * f_{xI, xc} * Ci[X_{I}] * 14)`"
   "S_I concentration", ":math:`S_{I, 1} = S_{I, in} + (f_{sI, xc} * X_{bio})`"
   "X_ch concentration", ":math:`X_{ch, 1} = X_{ch, in} + (f_{ch, xc} * X_{bio})`"
   "X_pr concentration", ":math:`X_{pr, 1} = X_{pr, in} + (f_{pr, xc} * X_{bio})`"
   "X_li concentration", ":math:`X_{li, 1} = X_{li, in} + (f_{li, xc} * X_{bio})`"
   "X_I concentration", ":math:`S_{I, 1} = X_{I, in} + (f_{xI, xc} * X_{bio})`"
   "S_IP concentration", ":math:`S_{IP, 1} = X_{PP, in} + (X_{bio} * Pi[X_{su}] * 31) - (X_{bio} * f_{sI, xc} * Pi[S_{I}] * 31) - (X_{bio} * f_{ch, xc} * P_ch) - (X_{bio} * f_{li, xc} * Pi[X_{li}] * 31) - (X_{bio} * f_{xI, xc} * Pi[X_{I}] * 31)`"
   "X_PHA concentration", ":math:`X_{PHA, 1} = 0`"
   "X_PP concentration", ":math:`X_{PP, 1} = 0`"
   "X_PAO concentration", ":math:`X_{PAO, 1} = 0`"
   "S_IC concentration", ":math:`S_{IC, 2} = (X_{ch, 1} * Ci[X_{ch}] * 12) + (X_{pr, 1} * Ci[X_{pr}] * 12) + (X_{li, 1} * Ci[X_{li}] * 12) - i_{CXS} * (X_{ch, 1} + X_{pr, 1} + X_{li, 1}) + (S_{su, in} * Ci[S_{su}] * 12) + (S_{aa, in} * Ci[S_{aa}] * 12) + (S_{fa, in} * Ci[S_{fa}] * 12) - i_{CSF} * (S_{su, in} + S_{aa, in} + S_{fa, in}) + (S_{va, in} * Ci[S_{va}] * 12) + (S_{bu, in} * Ci[S_{bu}] * 12) + (S_{pro, in} * Ci[S_{pro}] * 12) + (S_{AC, 1} * Ci[S_{ac}] * 12) - i_{CSA} * (S_{va, in} + S_{bu, in} + S_{pro, in} + S_{ac, 1})`"
   "S_IN concentration", ":math:`S_{IN, 2} = (X_{pr, 1} * Ni[X_{pr}] * 14)  - i_{NXS} * (X_{ch, 1} + X_{pr, 1} + X_{li, 1}) + (S_{aa, in} * Ni[S_{aa}] * 14) - i_{NSF} * (S_{su, in} + S_{fa, in} + S_{va, in})`"
   "S_IP concentration", ":math:`S_{IP, 2} = X_{PP, 1} + (X_{ch, 1} * P_ch) + (X_{li, 1} * Ni[X_{li}] * 14) - i_{PXS} * (X_{ch, 1} + X_{pr, 1} + X_{li, 1}) - i_{PSF} * (S_{su, in} + S_{fa, in} + S_{va, in})`"
   "S_K concentration", ":math:`S_{K, 2} = S_{K, in} + X_{PP, 1}`"
   "S_Mg concentration", ":math:`S_{Mg, 2} = S_{Mg, in} + X_{PP, 1}`"
   "S_F concentration", ":math:`S_{F, out} = S_{su, in} + S_{aa, in} + S_{fa, in}`"
   "S_A concentration", ":math:`S_{A, out} = S_{va, in} + S_{bu, in} + S_{pro, in} + S_{ac, 1}`"
   "S_I concentration", ":math:`S_{I, out} = S_{I, 1}`"
   "S_NH4 concentration", ":math:`S_{NH4, out} = S_{IN, in} + S_{IN, 1} + S_{IN, 2}`"
   "S_PO4 concentration", ":math:`S_{PO4, out} = S_{IP, in} + S_{IP, 1} + S_{IP, 2}`"
   "S_IC concentration", ":math:`S_{IC, out} = S_{IC, in} + S_{IC, 1} + S_{IC, 2}`"
   "X_I concentration", ":math:`X_{I, out} = X_{I, 1}`"
   "X_S concentration", ":math:`X_{S, out} = X_{ch, 1} + X_{pr, 1} + X_{li, 1}`"
   "S_K concentration", ":math:`S_{K, out} = S_{K, 2}`"
   "S_Mg concentration", ":math:`S_{Mg, out} = S_{Mg, 2}`"

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

