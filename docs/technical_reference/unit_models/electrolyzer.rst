Electrolyzer
===============================
This is a simplified electrolyzer unit model used to simulate electrolysis performance. The development of this unit model is ongoing.
   * simulation of this unit model is only supported with the Multi-Component Aqueous Solution (MCAS) property package
   * supports liquid phase only, vapor components are modeled in the liquid phase
   * supports steady-state only
   * assumes isothermal conditions, and performance is not temperature dependent
   * does not determine equilibrium of electrolysis products in solution
   * does not consider undesired reactions

.. index::
   pair: watertap.unit_models.electrolyzer;electrolyzer

.. currentmodule:: watertap.unit_models.electrolyzer

Introduction
------------
The basis for the electrolyzer model is considering the Faradaic conversion of species with respect to the electrolysis reactions at the cathode and anode. Faradaic conversion considers equating electrical current to the amount of electrons available for electrochemical reaction by Faraday's law.

Degrees of Freedom
------------------
In the core calculations of the electrolyzer unit model, there are 4 degrees of freedom in addition to the inlet state variables (i.e., temperature, pressure, component flowrates) that should be fixed for the model to be fully specified. In typical design cases the following 4 variables are typically fixed:

   * current, :math:`I`
   * current density, :math:`J`
   * current efficiency, :math:`\eta_{current}`
   * voltage efficiency, :math:`\eta_{voltage}`

Additionally, the electrolysis reactions at the anode and cathode must be specified. By default, the following stoichiometric variables will fix to 0. Therefore, only the nonzero stoichiometry must be specified. As an example for the chlor alkali process:

Anode:

    .. math::

        & Cl^- \to {\small\frac{1}{2}} Cl_2 + e^- \\\\

Cathode:

    .. math::

        & H_2O + e^- \to {\small\frac{1}{2}} H_2 + OH^- \\\\

The following variables should then be fixed:

.. code-block::

   # thermodynamic minimum voltage
   m.fs.electrolyzer.voltage_min.fix(2.2)

   # anode reactions
   m.fs.electrolyzer.anode_stoich["Liq", "CL-"].fix(-1)
   m.fs.electrolyzer.anode_stoich["Liq", "CL2-v"].fix(0.5)

   # cathode reactions
   m.fs.electrolyzer.cathode_stoich["Liq", "H2O"].fix(-1)
   m.fs.electrolyzer.cathode_stoich["Liq", "H2-v"].fix(0.5)
   m.fs.electrolyzer.cathode_stoich["Liq", "OH-"].fix(1)

   # membrane permeation (anode to cathode) satisfies charge balance
   m.fs.electrolyzer.membrane_stoich["Liq", "NA+"].fix(1)



Model Structure
------------------
The electrolyzer model consists of 2 ControlVolume0DBlocks, one for the anolyte and another for the catholyte. Currently, the generation of species via electrolysis reactions is handled by the `custom_molar_term` within the ControlVolume0DBlock.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "time", ":math:`t`", "[0]"
   "phases", ":math:`p`", "['Liq']"
   "components", ":math:`j`", "['H2O', solutes]"

.. _electrolyzer_variables:

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "molar generation of species j by electrolysis at the anode", ":math:`\dot{n}_{j,anode}`", "custom_reaction_anode", "t, j", ":math:`\frac{mol}{s}`"
   "molar generation of species j by electrolysis at the cathode", ":math:`\dot{n}_{j,cathode}`", "custom_reaction_cathode", "t, j", ":math:`\frac{mol}{s}`"
   "moles of electrons contributing to electrolysis", ":math:`\dot{n}_{e^-}`", "electron_flow", "None", ":math:`\frac{mol}{s}`"
   "current", ":math:`I`", "current", "None", ":math:`A`"
   "current efficiency", ":math:`\eta_{current}`", "efficiency_current", "None", ":math:`\text{dimensionless}`"
   "current density", ":math:`J`", "current_density", "None", ":math:`\frac{A}{m^3}`"
   "membrane area", ":math:`A_{membrane}`", "membrane_area", "None", ":math:`m^3`"
   "anode area", ":math:`A_{anode}`", "anode_area", "None", ":math:`m^3`"
   "cathode area", ":math:`A_{cathode}`", "cathode_area", "None", ":math:`m^3`"
   "applied voltage", ":math:`V`", "voltage_applied", "None", ":math:`V`"
   "voltage efficiency", ":math:`\eta_{voltage}`", "efficiency_voltage", "None", ":math:`V`"
   "power", ":math:`P`", "power", "None", ":math:`W`"
   "power efficiency", ":math:`\eta_{power}`", "efficiency_power", "None", ":math:`W`"

The following variables are constructed on the unit model for the current build. However, these variables are specific to the electrolysis reactions targeted.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "thermodynamic minimum voltage for electrolysis reactions to proceed", ":math:`V_{min}`", "voltage_min", "None", ":math:`V`"
   "stoichiometry of the reaction at the anode normalized to 1 electron", ":math:`\varepsilon_{j,anode}`", "anode_stoich", "p, j", ":math:`\text{dimensionless}`"
   "stoichiometry of the reaction at the cathode normalized to 1 electron", ":math:`\varepsilon_{j,cathode}`", "cathode_stoich", "p, j", ":math:`\text{dimensionless}`"
   "stoichiometry of the mass transfer across the membrane w.r.t. the electrolysis reaction", ":math:`\varepsilon_{j,membrane}`", "membrane_stoich", "p, j", ":math:`\text{dimensionless}`"

.. _electrolyzer_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "voltage efficiency as a function of minimum required by reaction", ":math:`V_{min} = V\eta_{voltage}`"
   "electrons passed between anode and cathode contributing to reactions", ":math:`\dot{n}_{e^-} = \frac{I\eta_{current}}{F}`"
   "membrane area", ":math:`I = JA_{membrane}`"
   "anode area", ":math:`I = JA_{anode}`"
   "cathode area", ":math:`I = JA_{cathode}`"
   "power", ":math:`P = IV`"
   "power efficiency", ":math:`P = \eta_{current}\eta_{voltage}`"
   "ion permeation through the membrane", ":math:`\dot{n}_j,membrane = -\varepsilon_{j,membrane}\dot{n}_{e^-}`"
   "molar generation of species according the anode electrolysis reaction", ":math:`\dot{n}_{j,anode} = \varepsilon_{j,anode}\dot{n}_{e^-}`"
   "molar generation of species according the cathode electrolysis reaction", ":math:`\dot{n}_{j,cathode} = \varepsilon_{j,cathode}\dot{n}_{e^-}`"

Costing Method
---------------

Costing Method Variables
+++++++++++++++++++++++++

The following parameters are constructed when applying the electrolyzer costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Default Value", "Units"

   "membrane replacement factor (fraction of membrane replaced/year)", ":math:`f_{membrane_replace}`", "factor_membrane_replacement", "0.33", ":math:`\text{dimensionless}`"
   "membrane unit cost", ":math:`c_{membrane}`", "membrane_unit_cost", "25", ":math:`\frac{$}{m^2}`"
   "anode unit cost", ":math:`c_{anode}`", "anode_unit_cost", "300", ":math:`\frac{$}{m^2}`"
   "cathode unit cost", ":math:`c_{cathode}`", "cathode_unit_cost", "600", ":math:`\frac{$}{m^2}`"
   "membrane, anode, and cathode fraction of total capital", ":math:`f_{material}`", "fraction_material_cost", "0.65", ":math:`\text{dimensionless}`"

The following variables are constructed when applying the electrolyzer costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units"

   "membrane capital cost", ":math:`C_{membrane}`", "membrane_cost", ":math:`$`"
   "anode capital cost", ":math:`C_{anode}`", "anode_cost", ":math:`$`"
   "cathode capital cost", ":math:`C_{cathode}`", "cathode_cost", ":math:`$`"
   "membrane replacement cost", ":math:`C_{membrane_replace}`", "membrane_replacement_cost", ":math:`frac{$}{yr}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital costs are contributing to the majority of material costs for the anode cathode and membrane. Each material cost is calculated individually then summed.

    .. math::

        & C_{membrane} = c_{membrane}A_{membrane} \\\\
        & C_{anode} = c_{anode}A_{anode} \\\\
        & C_{cathode} = c_{cathode}A_{cathode} \\\\
        & C_{cap,total} = \frac{C_{membrane}+C_{anode}+C_{cathode}}{f_{membrane_replace}}

Operating Cost Calculations
+++++++++++++++++++++++++++

Operating costs for the electrolyzer are the electricity requirements and membrane replacement costs. Electricity is costed using `cost_flow` applied to the `power` variable on the unit model.

    .. math::

        & C_{op,tot} = C_{membrane_replace} = f_{membrane_replace}c_{membrane}A_{membrane}

Code Documentation
-------------------

* :mod:`watertap.unit_models.electrolyzer`
* :mod:`watertap.costing.units.electrolyzer`

References
-----------
Hand, D. W., Crittenden, J. C., & Thacker, W. E. (1984). Simplified models for design of fixed-bed adsorption systems.
Journal of Environmental Engineering, 110(2), 440-456.

Crittenden, J., Rhodes, R., Hand, D., Howe, K., & Tchobanoglous, G. (2012). MWHs Water Treatment. Principles and Design.
John Wiley & Sons.

Crittenden, J. C., Berrigan, J. K., Hand, D. W., & Lykins, B. (1987). Design of Rapid Fixed‐Bed Adsorption Tests for
Nonconstant Diffusivities. Journal of Environmental Engineering, 113(2), 243–259.

United States Environmental Protection Agency. (2021). Work Breakdown Structure-Based Cost Model for Granular Activated
Carbon Drinking Water Treatment.