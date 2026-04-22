.. _electrolyzer:

Electrolyzer
============

.. code-block:: python

   from watertap.unit_models.electrolyzer import Electrolyzer

This is a simplified model used to approximate electrolyzer performance. With the current build, the model is simulated under the following assumptions:
   * simulation of this unit model is only supported with the :mod:`Multi-Component Aqueous Solution (MCAS) <watertap.property_models.multicomp_aq_sol_prop_pack>` property package
   * supports steady-state only
   * supports liquid phase only, vapor-liquid phase equilibrium is not calculated
   * does not determine chemical equilibrium of electrolysis products
   * assumes isothermal conditions and performance is not temperature-dependent
   * assumes isobaric conditions and performance is not pressure-dependent
   * does not consider the conversion of undesired reactions
   * does not incorporate the Nernst equation to calculate the potential as a function of the standard potential, the input of the electrochemical potential for the half-cell reactions is assumed to be provided at the operating conditions
   * conversion of reactants is predicted by a yield model, therefore rate laws and limitations of reaction and mass transfer are not evaluated


.. index::
   pair: watertap.unit_models.electrolyzer;electrolyzer

.. currentmodule:: watertap.unit_models.electrolyzer

Introduction
------------
This model was developed for the simulation of a chlor-alkali membrane electrolyzer, a conventional electrolyzer for the production of chlorine gas and hydroxide (Bommaraju & O’Brien, 2015) (Kent, 2007). This model has been demonstrated for the chlor-alkali configuration in the electrolyzer testing file. Though, it is structured such that it may be generalizable to other membrane electrolysis mechanisms, but it has not been validated for other processes.

With a membrane electrolyzer, the anolyte and catholyte are separated by an ion exchange membrane. The electrolyzers objective is dependent on the materials for the anode, cathode, and membrane. The anode and cathode are chosen to promote desired electrolysis reactions at the surface of the material. The membrane is chosen to selectively permeate ions between the electrolytes. Electrolyzer performance is largely quantified by the energy consumption, conversion of desired electrolysis products, and quality of the exiting anolyte and catholyte (Kumar et al., 2021). Faraday's law of electrolysis governs the relationship between electrical current supplied to the electrolytic cell and moles of electrons contributing to desired electrolysis reactions. The equation for Faradaic conversion is later established in "electrons contributing to reactions by Faraday's law of electrolysis" found in :ref:`Equations <electrolyzer_equations>`.

In this simplified model, the extent that the electrolysis reactions proceed is determined with respect to the ideal Faradaic conversion, resulting in a yield reaction model. As such, rate laws and limitations of reaction and mass transfer are not evaluated. For the efficiency of the reactions, the model defines a current efficiency variable (:math:`\theta_{current}`). Current efficiency describes the additional current required above the amount dictated by Faraday's law of electrolysis that contribute electrons to the desired electrolysis reaction. These hindrances are primarily the electrons that are supplied to the electrolyzer but contribute to undesired reactions. Since side reactions are not modeled in this simplified model, these electrons may be treated as unreacted. Additionally, this may also include generic loss of current (shunt currents) or neutralization of desired electrolysis products (back-migration across the ion exchange membrane or other) (Kumar et al., 2021). Nevertheless, the current efficiency describes electrons that are provided by a power supply, but do not contribute to the desired reaction. The model is structured that current efficiency is traditionally supplied as a fixed variable.

The second governing efficiency is the voltage efficiency (:math:`\theta_{voltage}`). Voltage efficiency describes the additional voltage above the thermodynamic minimum voltage (:math:`V_{rev}`) for the reaction to proceed (Kumar et al., 2021). This model directly defines these inefficiencies within the cell voltage equation. The voltage efficiency is then the contribution of anode overpotential, cathode overpotential, and ohmic resistance (related to voltage drop by Ohm's law) (Phillips, 2017). As the overpotentials and resistance are commonly known, these components of additional voltage are traditionally provided as the fixed variables. In this case, voltage efficiency is then calculated by constraints as a result of these components. However, if the intricacies of the electrolysis process are lesser known, this variable may be directly fixed and the model reworked to function given a fixed voltage efficiency.

Given that the power supplied is the product of voltage and current, a power efficiency may be defined that is the product of voltage efficiency and current efficiency. This power efficiency is then the ratio of the theoretical minimum power, under no inefficiencies, to the actual power supplied. In the unit model, the theoretical minimum power is not explicitly calculated. Where current efficiency and voltage efficiency are decoupled efficiencies that describe independent hindrances, power efficiency is the combined effect that directly related to energy demand and cost.

In this yield model, the surface area for electrochemical reaction on each electrode is the only relevant sizing variable required. This is used to validate that the current densities of each are tolerable by the material, and it is typically used as the fixed variable where the corresponding area is calculated as a degree of freedom. With no other sizing variables, effects in increased resistance through the electrolytes with further spaced electrodes and similar size dependent phenomena are not considered in the model.

Degrees of Freedom
------------------
For a given electrolyzer design, there are 8 degrees of freedom in addition to the inlet state variables (i.e., temperature, pressure, component flowrates) that should be fixed for the model to be fully specified. In typical design cases the following 8 variables are fixed:

   * membrane current density, :math:`i_{mem}`
   * anode current density, :math:`i_{ano}`
   * anode overpotential, :math:`\eta_{ano}`
   * cathode current density, :math:`i_{cat}`
   * cathode overpotential, :math:`\eta_{cat}`
   * current, :math:`I`
   * current efficiency, :math:`\theta_{current}`
   * resistance, :math:`R`

Additionally, the electrolysis reactions at the anode and cathode must be specified for their stoichiometry and electrochemical potential (Phillips, 2017). By default, the following stoichiometric variables will be fixed to 0 when the electrolyzer model block is constructed. Therefore, only the nonzero stoichiometry must be additionally specified. As an example for the chlor-alkali process, the following reactions occur with sodium ions permeating from the anolyte to catholyte. The stoichiometry is normalized to 1 mole of electrons as intended in the model framework. The electrochemical potential shown is the standard potential for the chlor-alkali reactions. In this model, the Nernst equation for calculating the potential as a function of the standard potential is not performed. Therefore, the input of the electrochemical potential for the half-cell reactions should be provided at the operating conditions of the electrolyzer.

Anode:

    .. math::

        & Cl^-\to{\small\frac{1}{2}}Cl_2+e^-\hspace{4em}E=1.21V \\\\

Cathode:

    .. math::

        & H_2O+e^-\to{\small\frac{1}{2}}H_2+OH^-\hspace{4em}E=-0.99V \\\\

The following variables should then be fixed:

.. code-block::

   # membrane properties
   m.fs.electrolyzer.membrane_ion_transport_number["Liq", "NA+"].fix(1)

   # anode properties
   # Cl- --> 0.5 Cl2 + e-
   m.fs.electrolyzer.anode_electrochem_potential.fix(1.21)
   m.fs.electrolyzer.anode_stoich["Liq", "CL-"].fix(-1)
   m.fs.electrolyzer.anode_stoich["Liq", "CL2-v"].fix(0.5)

   # cathode properties
   # H20 + e- --> 0.5 H2 + OH-
   m.fs.electrolyzer.cathode_electrochem_potential.fix(-0.99)
   m.fs.electrolyzer.cathode_stoich["Liq", "H2O"].fix(-1)
   m.fs.electrolyzer.cathode_stoich["Liq", "H2-v"].fix(0.5)
   m.fs.electrolyzer.cathode_stoich["Liq", "OH-"].fix(1)

Here, the flux of ions across the ion exchange membrane flux is handled by defining the membrane transport number. The sum of membrane transport numbers should correspond to satisfying the charge balance in both the anolyte and catholyte. The direction of membrane flow is calculated from anode to cathode in the unit model. Therefore a positive membrane transport number denotes flow of the corresponding ion (cation or anion depending on the ion exchange membrane) from the anolyte to catholyte.

Model Structure
------------------
The electrolyzer model consists of 2 ``ControlVolume0DBlocks``, one for the anolyte and another for the catholyte. Currently, the generation of species via electrolysis reactions is handled by the ``custom_molar_term`` within the ``ControlVolume0DBlock``. Coupling this custom molar term with the Faradaic conversion, the material balance may be written by the following and calculated within the ``ControlVolume0DBlock``. Considering the reaction block and reaction package are omitted, no temperature dependence and rigorous energy balance are considered.

Anode:

    .. math::

        & \dot{n}_{in,j}-\dot{n}_{out,j}+\frac{\left(\varepsilon_{ano,j}-t_{mem,j}\right)I\theta_{current}}{F}=0 \\\\

Cathode:

    .. math::

        & \dot{n}_{in,j}-\dot{n}_{out,j}+\frac{\left(\varepsilon_{cat,j}+t_{mem,j}\right)I\theta_{current}}{F}=0 \\\\

Scaling of unit model variables should first be performed using ``calculate_scaling_factors()``. For the case that the model is poorly scaled, the ``current`` variable should be rescaled. Estimated scaling factors are propagated from the ``current`` scaling factor for other unit model variables which are dependent on process scale. An example of the methodology is the provided below.

.. code-block::

   import idaes.core.util.scaling as iscale

   # fix the current for given design
   m.fs.electrolyzer.current.fix(2e8)

   # adjust the current's scaling factor to the inverse of the variable's magnitude
   iscale.set_scaling_factor(m.fs.electrolyzer.current, 1e-8)

   # run calculate scaling factor utility on model that propagates estimated scaling factors for other variables
   iscale.calculate_scaling_factors(m)

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

   "membrane area", ":math:`A_{mem}`", "membrane_area", "None", ":math:`m^2`"
   "membrane current density", ":math:`i_{mem}`", "membrane_current_density", "None", ":math:`\frac{A}{m^2}`"
   "ion transport number of species passing through the membrane :math:`^{ab}`", ":math:`t_{mem,j}`", "membrane_ion_transport_number", "[p, j]", ":math:`\text{dimensionless}`"
   "anode area", ":math:`A_{ano}`", "anode_area", "None", ":math:`m^2`"
   "anode current density", ":math:`i_{ano}`", "anode_current_density", "None", ":math:`\frac{A}{m^2}`"
   "anode electrochemical potential :math:`^a`", ":math:`E_{ano}`", "anode_electrochem_potential", "None", ":math:`V`"
   "anode overpotential", ":math:`\eta_{ano}`", "anode_overpotential", "None", ":math:`V`"
   "anode stoichiometry of the reaction :math:`^{ab}`", ":math:`\varepsilon_{ano,j}`", "anode_stoich", "[p, j]", ":math:`\text{dimensionless}`"
   "cathode area", ":math:`A_{cat}`", "cathode_area", "None", ":math:`m^2`"
   "cathode current density", ":math:`i_{cat}`", "cathode_current_density", "None", ":math:`\frac{A}{m^2}`"
   "cathode electrochemical potential :math:`^{ab}`", ":math:`E_{cat}`", "cathode_electrochem_potential", "None", ":math:`V`"
   "cathode overpotential", ":math:`\eta_{cat}`", "cathode_overpotential", "None", ":math:`V`"
   "cathode stoichiometry of the reaction :math:`^a`", ":math:`\varepsilon_{cat,j}`", "cathode_stoich", "[p, j]", ":math:`\text{dimensionless}`"
   "current", ":math:`I`", "current", "None", ":math:`A`"
   "cell voltage", ":math:`V_{cell}`", "voltage_cell", "None", ":math:`V`"
   "ohmic resistance", ":math:`R`", "resistance", "None", ":math:`\Omega`"
   "power", ":math:`P`", "power", "None", ":math:`W`"
   "reversible voltage", ":math:`V_{rev}`", "voltage_reversible", "None", ":math:`V`"
   "electrons contributing to electrolysis reactions", ":math:`\dot{n}_{e^-}`", "electron_flow", "None", ":math:`\frac{mol}{s}`"
   "current efficiency", ":math:`\theta_{current}`", "efficiency_current", "None", ":math:`\text{dimensionless}`"
   "voltage efficiency", ":math:`\theta_{voltage}`", "efficiency_voltage", "None", ":math:`\text{dimensionless}`"
   "power efficiency", ":math:`\theta_{power}`", "efficiency_power", "None", ":math:`\text{dimensionless}`"
   "molar flow of species j across the membrane from anolyte to catholyte", ":math:`\dot{n}_{mem,j}`", "mass_transfer_term", "[t, p, j]", ":math:`\frac{mol}{s}`"
   "molar generation of species j by electrolysis at the anode", ":math:`\dot{n}_{ano,j}`", "custom_reaction_anode", "[t, j]", ":math:`\frac{mol}{s}`"
   "molar generation of species j by electrolysis at the cathode", ":math:`\dot{n}_{cat,j}`", "custom_reaction_cathode", "[t, j]", ":math:`\frac{mol}{s}`"

| :math:`^a` Variable intended to be move to a interchangeable component block, callable to the base electrolyzer model
| :math:`^b` Value is normalized to 1 electron in electrolysis stoichiometry


For non-fixed efficiencies, electrochemical potentials, and other relevant variables, custom constraints  may be constructed at the flowsheet level. These may allow for predictive performance as a function of temperature, concentration, overpotential, and other variables.

.. _electrolyzer_equations:

Equations
-----------
.. csv-table::
   :header: "Description", "Equation"

   "membrane current density", ":math:`I = i_{mem}A_{mem}`"
   "ion permeation through the membrane", ":math:`\dot{n}_{mem,j} = -t_{mem,j}\dot{n}_{e^-}`"
   "anode current density", ":math:`I = i_{ano}A_{ano}`"
   "molar generation of species according the anode electrolysis reaction", ":math:`\dot{n}_{ano,j} = \varepsilon_{ano,j}\dot{n}_{e^-}`"
   "cathode current density", ":math:`I = i_{cat}A_{cat}`"
   "molar generation of species according the cathode electrolysis reaction", ":math:`\dot{n}_{cat,j} = \varepsilon_{cat,j}\dot{n}_{e^-}`"
   "reversible voltage", ":math:`V_{rev} = E_{ano}-E_{cat}`"
   "cell voltage", ":math:`V_{cell} = V_{rev}+\eta_{ano}+\eta_{cat}+IR`"
   "power", ":math:`P = IV_{cell}`"
   "electrons contributing to reactions by Faraday's law of electrolysis :math:`^c`", ":math:`\dot{n}_{e^-} = \frac{I\theta_{current}}{F}`"
   "voltage efficiency", ":math:`\theta_{voltage} = \frac{V_{rev}}{V_{cell}} = \frac{V_{rev}}{V_{rev}+\eta_{ano}+\eta_{cat}+IR}`"
   "power efficiency", ":math:`\theta_{power} = \theta_{current}\theta_{voltage}`"

\ :math:`^c` is Faraday's constant from ``idaes.core.util.constants``, 96,485 C/mol

References
-----------
Bommaraju, T. V., & O’Brien, T. F. (2015). Brine electrolysis. Electrochemistry Encyclopedia. https://knowledge.electrochem.org/encycl/art-b01-brine.htm

Kent, J. A. (Ed.). (2007). Kent and Riegel’s Handbook of Industrial Chemistry and Biotechnology. Springer US. https://doi.org/10.1007/978-0-387-27843-8

Kumar, A., Du, F., & Lienhard, J. H. (2021). Caustic Soda Production, Energy Efficiency, and Electrolyzers. ACS Energy Letters, 6(10), 3563–3566. https://doi.org/10.1021/acsenergylett.1c01827

O’Brien, T., Bommaraju, T. V., & Hine, F. (2005). Handbook of chlor-alkali technology. Springer.

Phillips, R., Edwards, A., Rome, B., Jones, D. R., & Dunnill, C. W. (2017). Minimising the ohmic resistance of an alkaline electrolysis cell through effective cell design. International Journal of Hydrogen Energy, 42(38), 23986–23994. https://doi.org/10.1016/j.ijhydene.2017.07.184

