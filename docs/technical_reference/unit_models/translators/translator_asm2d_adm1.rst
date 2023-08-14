ASM2d to ADM1 Translator
========================

Introduction
------------

A link is required to translate between biological based and physical or chemical mediated processes
to develop plant-wide modeling of wastewater treatment. This model mediates the interaction between
the Modified Activated Sludge Model 2d (ASM2d) and the Modified Anaerobic Digestor Model 1 (ADM1).

The model relies on the following key assumptions:

   * supports only liquid phase
   * supports only Modified ASM2d to Modified ADM1 translations

.. index::
   pair: watertap.unit_models.translators.translator_adm1_asm2d;translator_asm2d_adm1

.. currentmodule:: watertap.unit_models.translators.translator_asm2d_adm1

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
   "Inlet Components", ":math:`j_{in}`", "['H2O', 'S_A', 'S_F', 'S_I', 'S_N2', 'S_NH4', 'S_NO3', 'S_O2', 'S_PO4', 'S_K', 'S_Mg', 'S_IC', 'X_AUT', 'X_H', 'X_I', 'X_PAO', 'X_PHA', 'X_PP', 'X_S']"
   "Outlet Components", ":math:`j_{out}`", "['H2O', 'S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_IP', 'S_I', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'X_PHA', 'X_PP', 'X_PAO', 'S_K', 'S_Mg']"
   "Ion", ":math:`j_{in}`", "['S_cat', 'S_an'] \  :sup:`1`"
   "Zero Flow Components", ":math:`z`", "['S_fa', 'S_h2', 'S_ch4', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2']"

**Notes**
 :sup:`1` "Ion" is a subset of "Outlet Components" and uses the same symbol j_in.

Parameters
----------

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Value", "Units"

   "Soluble inerts from composites", ":math:`f_{sI, xc}`", "f_sI_xc", 1e-9, ":math:`\text{dimensionless}`"
   "Particulate inerts from composites", ":math:`f_{xI, xc}`", "f_xI_xc", 0.1, ":math:`\text{dimensionless}`"
   "Carbohydrates from composites", ":math:`f_{ch, xc}`", "f_ch_xc", 0.275, ":math:`\text{dimensionless}`"
   "Proteins from composites", ":math:`f_{pr, xc}`", "f_pr_xc", 0.275, ":math:`\text{dimensionless}`"
   "Lipids from composites", ":math:`f_{li, xc}`", "f_li_xc", 0.35, ":math:`\text{dimensionless}`"
   "Valerate from polyhydroxyalkanoates", ":math:`f_{XPHA, Sva}`", "f_XPHA_Sva", 0.1, ":math:`\text{dimensionless}`"
   "Butyrate from polyhydroxyalkanoates", ":math:`f_{XPHA, Sbu}`", "f_XPHA_Sbu", 0.1, ":math:`\text{dimensionless}`"
   "Propionate from polyhydroxyalkanoates", ":math:`f_{XPHA, Spro}`", "f_XPHA_Spro", 0.4, ":math:`\text{dimensionless}`"
   "Acetate from polyhydroxyalkanoates", ":math:`f_{XPHA, Sac}`", "f_XPHA_Sac", 0.4, ":math:`\text{dimensionless}`"
   "Carbon content of polyhydroxyalkanoates", ":math:`C_{PHA}`", "C_PHA", 0.025, ":math:`\text{dimensionless}`"


.. _Translator_ASM2d_ADM1_equations:

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Volumetric flow equality", ":math:`F_{out} = F_{in}`"
   "Temperature balance", ":math:`T_{out} = T_{in}`"
   "Pressure balance", ":math:`P_{out} = P_{in}`"
   "Zero-flow component conversions", ":math:`C_{z, out} = 0`"
   "Anions balance", ":math:`S_{an} = \frac{S_{IN, out}}{14}`"
   "Cations balance", ":math:`S_{cat} = \frac{S_{IC, out}}{12}`"
   "COD demanding compounds in S_O2", ":math:`COD_{SO2} = S_{O2, in} / \frac{1 - Y_{H}}{Y_{H}}`"
   "S_O2 concentration", ":math:`S_{O2, 1} = S_{O2, in} - \frac{1 - Y_{H}}{Y_{H}} * COD_{SO2}`"
   "S_A concentration", ":math:`S_{A, 1} = S_{A, in} - \frac{COD_{SO2}}{Y_{H}}`"
   "S_NH4 concentration", ":math:`S_{NH4, 1} = S_{NH4, in} - (i_{NBM} * COD_{SO2})`"
   "S_PO4 concentration", ":math:`S_{PO4, 1} = S_{PO4, in} - (i_{PBM} * COD_{SO2})`"
   "S_IC concentration", ":math:`S_{IC, 1} = S_{IC, in} + \frac{COD_{SO2} * i_{CSA}}{Y_{H}} + (COD_{SO2} * i_{CXB})`"
   "X_H concentration", ":math:`X_{H, 1} = X_{H, in} + COD_{SO2}`"
   "COD demanding compounds in S_NO3", ":math:`COD_{SNO3} = S_{NO3, in} / \frac{1 - Y_{H}}{i_{NOx, N2} * Y_{H}}`"
   "S_A concentration", ":math:`S_{A, 2} = S_{A, 1} - \frac{COD_{SNO3}}{Y_{H}}`"
   "S_NH4 concentration", ":math:`S_{NH4, 2} = S_{NH4, 1} - (COD_{SNO3} * i_{NBM})`"
   "S_N2 concentration", ":math:`S_{N2, 2} = S_{N2, in} + \frac{1 - Y_{H}}{i_{NOx, N2} * Y_{H}} * COD_{SNO3}`"
   "S_NO3 concentration", ":math:`S_{NO3, 2} = S_{NO3, in} - \frac{1 - Y_{H}}{i_{NOx, N2} * Y_{H}} * COD_{SNO3}`"
   "S_PO4 concentration", ":math:`S_{PO4, 2} = S_{PO4, 1} - (COD_{SNO3} * i_{PBM})`"
   "S_IC concentration", ":math:`S_{IC, 2} = S_{IC, 1} + \frac{COD_{SNO3} * i_{CSA}}{Y_{H}} + (COD_{SNO3} * i_{CXB})`"
   "X_H concentration", ":math:`X_{H, 2} = X_{H, 1} + COD_{SNO3}`"
   "Nitrogen demand for soluble inerts", ":math:`S_{ND} = S_{F, in} * i_{NSF}`"
   "Phosphorus demand for soluble inerts", ":math:`S_{PD} = S_{F, in} * i_{PSF}`"
   "Organic nitrogen from soluble inerts", ":math:`SN_{org} = \frac{S_{ND}}{Ni[S_{aa}] * 14}`"
   "Monosaccharides mapping (if :math:`SN_{org} >= S_{F, in}`)", ":math:`S_{su} = in`"
   "Monosaccharides mapping (if :math:`SN_{org} < S_{F, in}`)", ":math:`S_{su} = \frac{S_{F, in} - SN_{org}}{1000}`"
   "Amino acids mapping (if :math:`SN_{org} >= S_{F, in}`)", ":math:`S_{aa} = \frac{S_{F, in}}{1000}`"
   "Amino acids mapping (if :math:`SN_{org} < S_{F, in}`)", ":math:`S_{aa} = \frac{SN_{org}}{1000}`"
   "S_F concentration", ":math:`S_{F, 3} = S_{F, in} - (S_{su} * 1000) - (S_{aa} * 1000)`"
   "S_NH4 concentration", ":math:`S_{NH4, 3} = S_{NH4, 2} + (S_{F, in} * i_{NSF}) - (S_{aa} * Ni[S_{aa}] * 1000 * 14)`"
   "S_PO4 concentration", ":math:`S_{PO4, 3} = S_{PO4, 2} + (S_{F, in} * i_{PSF})`"
   "S_IC concentration", ":math:`S_{IC, 3} = S_{IC, 2} + (S_{F, in} * i_{CSF}) - (S_{su} * Ci[S_{su}] * 1000 * 12) - (S_{aa} * Ci[S_{aa}] * 1000 * 12)`"

Equations and Relationships (With Decay)
----------------------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Biomass concentration", ":math:`bio = X_{H, 2} + X_{PAO, in} + X_{AUT, in}`"
   "S_I concentration", ":math:`S_{I, 4} = S_{I, in} + (bio * f_{sI, xc})`"
   "S_NH4 concentration", ":math:`S_{NH4, 4} = S_{NH4, 3} + (bio * i_{NBM}) - (bio * f_{sI, xc} * i_{NSI}) - (bio * f_{xI, xc} * i_{NSI}) - (bio * f_{pr, xc} * Ni[X_{pr}] * 14)`"
   "S_PO4 concentration", ":math:`S_{PO4, 4} = S_{PO4, 3} + (bio * i_{PBM}) - (bio * f_{sI, xc} * i_{PSI}) - (bio * f_{xI, xc} * i_{PXI}) - (bio * f_{ch, xc} * Pi[X_{ch}] * 31) - (bio * f_{li, xc} * Pi[X_{li}] * 31)`"
   "S_IC concentration", ":math:`S_{IC, 4} = S_{IC, 3} + (bio * i_{CXB}) - (bio * f_{sI, xc} * i_{CSI}) - (bio * f_{xI, xc} * i_{CXI}) - (bio * f_{pr, xc} * Ci[X_{pr}] * 12) - (bio * f_{ch, xc} * Ci[X_{ch}] * 12) - (bio * f_{li, xc} * Ci[X_{li}] * 12)`"
   "X_I concentration", ":math:`X_{I, 4} = X_{I, in} + (bio * f_{xI, xc})`"
   "X_H concentration", ":math:`X_{H, 4} = 0`"
   "X_PAO concentration", ":math:`X_{PAO, 4} = 0`"
   "X_AUT concentration", ":math:`X_{AUT, 4} = 0`"
   "Nitrogen demand for particulate inerts", ":math:`X_{ND} = X_{S, in} * i_{NXS}`"
   "Phosphorus demand for particulate inerts", ":math:`X_{PD} = X_{S, in} * i_{PXS}`"
   "Organic nitrogen from particulate inerts", ":math:`XN_{org} = \frac{X_{ND}}{Ni[X_{pr}] * 14}`"
   "Carbohydrates mapping (if :math:`XN_{org} >= X_{S, in}`)", ":math:`X_{ch} = 0`"
   "Carbohydrates mapping (if :math:`XN_{org} < X_{S, in}`)", ":math:`X_{ch} = \frac{(X_{S, in} - XN_{org}) * 0.4}{1000}`"
   "Protein mapping (if :math:`XN_{org} >= X_{S, in}`)", ":math:`X_{pr} = \frac{S_{F, 3}}{1000}`"
   "Protein mapping (if :math:`XN_{org} < X_{S, in}`)", ":math:`X_{pr} = \frac{XN_{org}}{1000}`"
   "Lipids mapping (if :math:`XN_{org} >= X_{S, in}`)", ":math:`X_{li} = 0`"
   "Lipids mapping (if :math:`XN_{org} < X_{S, in}`)", ":math:`X_{li} = \frac{(X_{S, in} - XN_{org}) * 0.6}{1000}`"
   "S_NH4 concentration", ":math:`S_{NH4, 5} = S_{NH4, 4} + (X_{S, in} * i_{NXS}) - (X_{pr} * Ni[X_{pr}] * 1000 * 14)`"
   "S_PO4 concentration", ":math:`S_{PO4, 5} = S_{PO4, 4} + (X_{S, in} * i_{PXS}) - (X_{ch} * Pi[X_{ch}] * 1000 * 31) - (X_{li} * Pi[X_{li}] * 1000 * 31)`"
   "S_IC concentration", ":math:`S_{IC, 5} = S_{IC, 4} + (S_{F, in} * i_{CXS}) - (X_{ch} * Ci[X_{ch}] * 1000 * 12) - (X_{pr} * Ci[X_{pr}] * 1000 * 12) - (X_{li} * Ci[X_{li}] * 1000 * 12)`"
   "X_S concentration", ":math:`X_{S, 5} = 0`"
   "X_PP concentration", ":math:`X_{PP, 6} = 0`"
   "X_PHA concentration", ":math:`X_{PHA, 6} = 0`"
   "S_va concentration", ":math:`S_{va, 6} = \frac{X_{PHA, 6} * f_{XPHA, Sva}}{1000}`"
   "S_bu concentration", ":math:`S_{bu, 6} = \frac{X_{PHA, 6} * f_{XPHA, Sbu}}{1000}`"
   "S_pro concentration", ":math:`S_{pro, 6} = \frac{X_{PHA, 6} * f_{XPHA, Spro}}{1000}`"
   "S_ac concentration", ":math:`S_{ac, 6} = \frac{X_{PHA, 6} * f_{XPHA, Sac}}{1000}`"
   "S_PO4 concentration", ":math:`S_{PO4, 6} = S_{PO4, 5} + X_{PP, in}`"
   "S_IC concentration", ":math:`S_{IC, 6} = S_{IC, 5} + (X_{PHA, in} * C_{PHA}) - (S_{va, 6} * Ci[S_{va}] * 1000 * 12) - (S_{bu, 6} * Ci[S_{bu}] * 1000 * 12) - (S_{pro, 6} * Ci[S_{pro}] * 1000 * 12) - (S_{ac, 6} * Ci[S_{ac}] * 1000 * 12)`"
   "S_K concentration", ":math:`S_{K, 6} = S_{K, in} + (K_{XPP} * X_{PP, in})`"
   "S_Mg concentration", ":math:`S_{Mg, 6} = S_{Mg, in} + (Mg_{XPP} * X_{PP, in})`"

Equations and Relationships (Without Decay)
-------------------------------------------

.. csv-table::
   :header: "Description", "Equation"

   "Biomass concentration", ":math:`bio = X_{H, 2} + X_{AUT, in}`"
   "S_I concentration", ":math:`S_{I, 4} = S_{I, in} + (bio * f_{sI, xc})`"
   "S_NH4 concentration", ":math:`S_{NH4, 4} = S_{NH4, 3} + (bio * i_{NBM}) - (bio * f_{sI, xc} * i_{NSI}) - (bio * f_{xI, xc} * i_{NSI}) - (bio * f_{pr, xc} * Ni[X_{pr}] * 14)`"
   "S_PO4 concentration", ":math:`S_{PO4, 4} = S_{PO4, 3} + (bio * i_{PBM}) - (bio * f_{sI, xc} * i_{PSI}) - (bio * f_{xI, xc} * i_{PXI}) - (bio * f_{ch, xc} * Pi[X_{ch}] * 31) - (bio * f_{li, xc} * Pi[X_{li}] * 31)`"
   "S_IC concentration", ":math:`S_{IC, 4} = S_{IC, 3} + (bio * i_{CXB}) - (bio * f_{sI, xc} * i_{CSI}) - (bio * f_{xI, xc} * i_{CXI}) - (bio * f_{pr, xc} * Ci[X_{pr}] * 12) - (bio * f_{ch, xc} * Ci[X_{ch}] * 12) - (bio * f_{li, xc} * Ci[X_{li}] * 12)`"
   "X_I concentration", ":math:`X_{I, 4} = X_{I, in} + (bio * f_{xI, xc})`"
   "X_H concentration", ":math:`X_{H, 4} = 0`"
   "X_PAO concentration", ":math:`X_{PAO, 4} = X_{PAO, in}`"
   "X_PP concentration", ":math:`X_{PP, 4} = X_{PP, in}`"
   "X_PHA concentration", ":math:`X_{PHA, 4} = X_{PHA, in}`"
   "X_AUT concentration", ":math:`X_{AUT, 4} = 0`"
   "Nitrogen demand for particulate inerts", ":math:`X_{ND} = X_{S, in} * i_{NXS}`"
   "Phosphorus demand for particulate inerts", ":math:`X_{PD} = X_{S, in} * i_{PXS}`"
   "Organic nitrogen from particulate inerts", ":math:`XN_{org} = \frac{X_{ND}}{Ni[X_{pr}] * 14}`"
   "Carbohydrates mapping (if :math:`XN_{org} >= X_{S, in}`)", ":math:`X_{ch} = 0`"
   "Carbohydrates mapping (if :math:`XN_{org} < X_{S, in}`)", ":math:`X_{ch} = \frac{(X_{S, in} - XN_{org}) * 0.4}{1000}`"
   "Protein mapping (if :math:`XN_{org} >= X_{S, in}`)", ":math:`X_{pr} = \frac{S_{F, 3}}{1000}`"
   "Protein mapping (if :math:`XN_{org} < X_{S, in}`)", ":math:`X_{pr} = \frac{XN_{org}}{1000}`"
   "Lipids mapping (if :math:`XN_{org} >= X_{S, in}`)", ":math:`X_{li} = 0`"
   "Lipids mapping (if :math:`XN_{org} < X_{S, in}`)", ":math:`X_{li} = \frac{(X_{S, in} - XN_{org}) * 0.6}{1000}`"
   "S_NH4 concentration", ":math:`S_{NH4, 5} = S_{NH4, 4} + (X_{S, in} * i_{NXS}) - (X_{pr} * Ni[X_{pr}] * 1000 * 14)`"
   "S_PO4 concentration", ":math:`S_{PO4, 5} = S_{PO4, 4} + (X_{S, in} * i_{PXS}) - (X_{ch} * Pi[X_{ch}] * 1000 * 31) - (X_{li} * Pi[X_{li}] * 1000 * 31)`"
   "S_IC concentration", ":math:`S_{IC, 5} = S_{IC, 4} + (S_{F, in} * i_{CXS}) - (X_{ch} * Ci[X_{ch}] * 1000 * 12) - (X_{pr} * Ci[X_{pr}] * 1000 * 12) - (X_{li} * Ci[X_{li}] * 1000 * 12)`"
   "X_S concentration", ":math:`X_{S, 5} = 0`"
   "X_PAO concentration", ":math:`X_{PAO, 5} = X_{PAO, in}`"
   "X_PP concentration", ":math:`X_{PP, 5} = X_{PP, in}`"
   "X_PHA concentration", ":math:`X_{PHA, 5} = X_{PHA, in}`"
   "S_va concentration", ":math:`S_{va, 6} = 0`"
   "S_bu concentration", ":math:`S_{bu, 6} = 0`"
   "S_pro concentration", ":math:`S_{pro, 6} = 0`"
   "S_ac concentration", ":math:`S_{ac, 6} = 0`"
   "X_PAO concentration", ":math:`X_{PAO, 6} = X_{PAO, in}`"
   "X_PP concentration", ":math:`X_{PP, 6} = X_{PP, in}`"
   "X_PHA concentration", ":math:`X_{PHA, 6} = X_{PHA, in}`"
   "S_va concentration", ":math:`S_{va, 6} = \frac{X_{PHA, 6} * f_{XPHA, Sva}}{1000}`"
   "S_bu concentration", ":math:`S_{bu, 6} = \frac{X_{PHA, 6} * f_{XPHA, Sbu}}{1000}`"
   "S_pro concentration", ":math:`S_{pro, 6} = \frac{X_{PHA, 6} * f_{XPHA, Spro}}{1000}`"
   "S_ac concentration", ":math:`S_{ac, 6} = \frac{X_{PHA, 6} * f_{XPHA, Sac}}{1000}`"
   "S_PO4 concentration", ":math:`S_{PO4, 6} = S_{PO4, 5} + X_{PP, in}`"
   "S_IC concentration", ":math:`S_{IC, 6} = S_{IC, 5} + (X_{PHA, in} * C_{PHA}) - (S_{va, 6} * Ci[S_{va}] * 1000 * 12) - (S_{bu, 6} * Ci[S_{bu}] * 1000 * 12) - (S_{pro, 6} * Ci[S_{pro}] * 1000 * 12) - (S_{ac, 6} * Ci[S_{ac}] * 1000 * 12)`"
   "S_K concentration", ":math:`S_{K, 6} = S_{K, in} + (K_{XPP} * X_{PP, in})`"
   "S_Mg concentration", ":math:`S_{Mg, 6} = S_{Mg, in} + (Mg_{XPP} * X_{PP, in})`"


Classes
-------
.. currentmodule:: watertap.unit_models.translators.translator_asm2d_adm1

.. autoclass:: TranslatorDataASM2dADM1
    :members:
    :noindex:

References
----------
[1] Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J., 2016.
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes.
Water Research, 95, pp.370-382. https://github.com/wwtmodels/Plant-Wide-Models

