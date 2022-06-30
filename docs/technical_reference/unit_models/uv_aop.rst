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
   "Mass flowrate of components", ":math:`M_j`", "flow_mass_phase_comp", "[t, 'Liq', j]", ":math:`\text{kg/s}`"
   "Inlet volumetric flowrate", ":math:`F_{in}`", "flow_vol", "[t]", ":math:`\text{m}^3\text{/s}`"

**NOTE: Variables for 'temperature', 'pressure', 'flow_mass_phase_comp', and 'flow_vol' come from the associated property package as state variables and are accessed via {port_name}.{state_var_name}**

Aside from the inlet feed state variables (i.e., temperature, pressure, component mass flowrates),
the UV AOP model has at least an additional 5 degrees of freedom that
the user must specify. The table below gives an outline of these.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Inactivation rate coefficient", ":math:`k`", "inactivation_rate", "['Liq', j]", ":math:`\text{m}^2\text{/J}`"
   "Pseudo-first order rate constant", ":math:`k_0`", "rate_constant", "['Liq', j]", ":math:`\text{s}^{-1}`"
   "UV dose", ":math:`D`", "uv_dose", None, ":math:`\text{J/}\text{m}^2`"
   "Average intensity of UV light", ":math:`I`", "uv_intensity", None, ":math:`\text{J/}\text{m}^2\text{/s}`"
   "Exposure time of UV light", ":math:`t`", "exposure_time", None, ":math:`\text{s}`"
   "Electricity demand of components", ":math:`E_j`", "electricity_demand_phase_comp", "[t, 'Liq', j]", ":math:`\text{W}`"
   "Electricity efficiency per log order reduction (EE/O)", ":math:`EE/O_j`", "electrical_efficiency", "[t, 'Liq', j]", ":math:`\text{J/}\text{m}^3`"
   "Lamp efficiency", ":math:`\eta`", "lamp_efficiency", None, None

**Users must provide values for and 'fix' certain variables to solve the model with DOF=0. Thus, users should fix
    * either 'inactivation_rate' and 'rate_constant',

    * two variables out of 'uv_dose', 'uv_intensity' and 'exposure_time',
    * either 'electricity_demand_phase_comp' and 'electrical_efficiency', 

    * and 'lamp_efficiency'.


However, users may later unfix certain variables for optimization purposes.**

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "UV dose", ":math:`D = I \cdot t`"
   "Pseudo-first order rate constant", ":math:`k_0 = I \cdot k`"
   "Solvent mass balance", ":math:`M_{\text{H2O},out} = M_{\text{H2O},in}`"
   "Solute mass balance", ":math:`M_{j,out} = M_{j,in} \cdot \exp(D \cdot k)`"
   "Electricity demand", ":math:`E_j = EE/O_j \cdot F_{in} \cdot \log_{10}(M_{j,in} / M_{j,out}) / \eta`"

Class Documentation
-------------------

* :mod:`watertap.unit_models.uv_aop`

References
----------
Wright, H., Gaithuma, D., Greene, D., Aieta, M. (2006) Integration of validation, design, and operation provides optimal implementation of UV disinfection.
American Water Works Association, 98, 81-92