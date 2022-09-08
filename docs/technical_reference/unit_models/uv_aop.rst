Ultraviolet Advanced Oxidation Process
======================================

.. index::
   pair: watertap.unit_models.uv_aop;uv_aop

.. currentmodule:: watertap.unit_models.uv_aop

Introduction
------------

Ultraviolet (UV) light is widely used in the water industry to inactivate microorganisms by damaging their nucleic acid, which prevents the pathogen from replicating and causing infection (Wright et al, 2006). Some contaminants resistant to biodegradation, such as N-Nitrosodimethylamine (NDMA), are susceptible to UV disinfection. The high oxidizing capability of UV is often used prior to chemical disinfection to ensure high microbial quality water as well as significantly reducing the chemical dosage required.

Advanced oxidation processes (AOPs) are technologies involving the generation of highly reactive oxidative species, predominantly hydroxyl radicals (:math:`\text{HO} \cdot`). Unlike conventional chemical oxidation processes, such as using chlorine, which are selective as to which compounds they can degrade, AOPs are able to completely convert organic compounds into carbon dioxide, water and mineral acids. In addition, AOPs are feasible for full-scale use to destroy organic compounds because they generate hydroxyl radicals at ambient temperature and atmospheric pressure. Typically, the commercially available UV AOPs for industrial water treatment are:

1) UV light and hydrogen peroxide (:math:`\text{H}_2\text{O}_2`)

2) UV light and ozone (:math:`\text{O}_3`)

3) UV light, hydrogen peroxide (:math:`\text{H}_2\text{O}_2`) and ozone (:math:`\text{O}_3`)

UV AOPs can be  modeled at several different levels, depending on the known kinetic pathways and the modeling objectives. In this work, a basic level kinetic model is presented with an assumption on pseudo-steady state approximation for the kinetic description of free radical species. A pseudo-first order rate constant is utilized to represent the overall degradation rate of contaminants. The users need to provide either UV dose and disinfection rate, or rate constant and exposure time that are can be acquired during UV validation tests. These measurements are then used to simulate the removal of contaminants during the UV AOP process.

This model also accounts for the energy demand for the disinfection process. The users need to provide Electrical Efficiency per Log Order Reduction (EE/O) to evaluate energy costs.

Ports
-----

The model provides two ports (Pyomo notation in parentheses):

* Inlet port (inlet)
* Outlet port (outlet)

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NDMA', ...]"

Degrees of Freedom and Variables
--------------------------------
The UV system includes the state variables from the associated property package, and an outline of these variables are listed below:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inlet temperature", ":math:`T`", "temperature", "[t]", ":math:`\text{K}`"
   "Inlet pressure", ":math:`P_{in}`", "pressure", "[t]", ":math:`\text{Pa}`"
   "Outlet pressure", ":math:`P_{out}`", "pressure", "[t]", ":math:`\text{Pa}`"
   "Mass flowrate of components", ":math:`M_{p, j}`", "flow_mass_phase_comp", "[t, p, j]", ":math:`\text{kg/s}`"
   "Inlet volumetric flowrate", ":math:`F_{in}`", "flow_vol", "[t]", ":math:`\text{m}^3\text{/s}`"

**NOTE: Variables for 'temperature', 'pressure', 'flow_mass_phase_comp', and 'flow_vol' come from the associated property package as state variables and are accessed via {port_name}.{state_var_name}**

Aside from the inlet feed state variables (i.e., temperature, pressure, component mass flowrates),
the UV AOP model has at least an additional 6 degrees of freedom that
the user must specify when the `uv_dose_type` configuration option is set to `UVDoseType.fixed`.
The table below gives an outline of these.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inactivation rate coefficient", ":math:`k`", "inactivation_rate", "[p, j]", ":math:`\text{m}^2\text{/J}`"
   "Overall pseudo-first order rate constant", ":math:`k_0`", "rate_constant", "[p, j]", ":math:`\text{s}^{-1}`"
   "Pseudo-first order rate constant for direct photolysis of component", ":math:`k_d`", "photolysis_rate_constant", "[p, j]", ":math:`\text{s}^{-1}`"
   "Pseudo-first order rate constant for indirect photolysis of component", ":math:`k_i`", "reaction_rate_constant", "[p, j]", ":math:`\text{s}^{-1}`"
   "UV dose", ":math:`D`", "uv_dose", None, ":math:`\text{J/}\text{m}^2`"
   "Average intensity of UV light", ":math:`I`", "uv_intensity", None, ":math:`\text{J/}\text{m}^2\text{/s}`"
   "Exposure time of UV light", ":math:`t`", "exposure_time", None, ":math:`\text{s}`"
   "Electricity demand of components with phase", ":math:`E_{p, j}`", "electricity_demand_phase_comp", "[t, p, j]", ":math:`\text{W}`"
   "Electricity demand of components", ":math:`E_j`", "electricity_demand_comp", "[t, j]", ":math:`\text{W}`"
   "Electricity demand", ":math:`E`", "electricity_demand", "[t]", ":math:`\text{W}`"
   "Electricity efficiency per log order reduction (EE/O)", ":math:`EE/O_{p, j}`", "electrical_efficiency_phase_comp", "[t, p, j]", ":math:`\text{J/}\text{m}^3`"
   "Lamp efficiency", ":math:`\eta`", "lamp_efficiency", None, None

**Users must provide values for and 'fix' certain variables to solve the model with DOF=0. Thus, users should fix**
    * either 'inactivation_rate' or 'rate_constant',
    * either 'photolysis_rate_constant' or 'reaction_rate_constant',
    * two variables out of 'uv_dose', 'uv_intensity' and 'exposure_time',
    * either 'electricity_demand_phase_comp' or 'electrical_efficiency',
    * and 'lamp_efficiency'.

**However, users may later unfix certain variables for optimization purposes.**

When setting the `uv_dose_type` configuration option to `UVDoseType.calculated`, there are 7 additional variables that must be fixed. This leads to a total of 13 degrees of freedom. Additional variables that must be fixed include:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "UV transmittance", ":math:`UVT`", "UVT", None, None
   "Coefficient A", ":math:`A`", "A_coeff", None, None
   "Coefficient B", ":math:`B`", "B_coeff", None, None
   "Coefficient C", ":math:`C`", "C_coeff", None, None
   "Coefficient D", ":math:`D`", "D_coeff", None, None
   "Relative lamp output", ":math:`\frac{S}{S_0}`", "relative_lamp_output", None, None
   "Number of banks", ":math:`N_{bank}`", "num_of_banks", None, None

Equations and Relationships
---------------------------
if ``uv_dose_type`` is set to default ``UVDoseType.fixed``:

.. csv-table::
   :header: "Description", "Equation"

   "UV dose", ":math:`D = I \cdot t`"
   "Pseudo-first order rate constant", ":math:`k_0 = I \cdot k`"
   "Pseudo-first order rate constant", ":math:`k_0 = k_d + k_i`"
   "Solvent mass balance", ":math:`M_{\text{H2O},out} = M_{\text{H2O},in}`"
   "Solute mass balance", ":math:`M_{p, j, out} = M_{p, j, in} \cdot \exp(D \cdot k)`"
   "Electricity demand of each component with phase", ":math:`E_{p, j} = EE/O_{p, j} \cdot F_{in} \cdot \log_{10}(M_{p, j, in} / M_{p, j, out}) / \eta`"
   "Electricity demand of each component", ":math:`E_j = \max_p E_{p, j}`"
   "Electricity demand", ":math:`E = \max_j E_j`"

if ``uv_dose_type`` is set to ``UVDoseType.calculated``, there is one additional equation:

.. csv-table::
   :header: "Description", "Equation"

   "UV dose", ":math:`D = 10^A \cdot (-\log_{10}(UVT))^{(-B \cdot \log_{10}(UVT))} \cdot (\frac{S}{S_0} / F_{in})^C \cdot N_{bank}^D`"

Class Documentation
-------------------

* :mod:`watertap.unit_models.uv_aop`

References
----------
Wright, H., Gaithuma, D., Greene, D., Aieta, M. (2006) Integration of validation, design, and operation provides optimal implementation of UV disinfection.
American Water Works Association, 98, 81-92