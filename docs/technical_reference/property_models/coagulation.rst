Coagulation Property Package
============================

This package implements property relationships for water density as a function of
temperature, pressure, and mass fraction of suspended/dissolved solids from
`Engineering Toolbox. (2003) <https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html>`_
and water viscosity as a function of temperature from
`D.S. Viswananth, G. Natarajan. (1989) <https://www.osti.gov/biblio/6562161>`_.

**Note: TDS = Total Dissolved Solids and TSS = Total Suspended Solids**

This coagulation property package:
   * supports only 'H2O', 'TDS', 'TSS', and 'Sludge' as Components
   * supports only liquid phase
   * is formulated on a mass basis
   * does NOT support formulations on a molar basis
   * includes mass density correction for fraction of suspended/dissolved solids

Sets
----
.. csv-table::
  :header: "Description", "Symbol", "Indices"

  "Components", ":math:`j`", "['H2O', 'TDS', 'TSS', 'Sludge']"
  "Phases", ":math:`p`", "['Liq']"

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mass flowrate", ":math:`M_j`", "flow_mass_phase_comp", "[p, j]", ":math:`\text{kg/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"

Parameters
----------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Index", "Units"

 "Fluid specific heat capacity", ":math:`c_p`", "cp", "None", ":math:`\text{J/kg/K}`"
 "Reference density (at 273 K)", ":math:`\rho_{ref}`", "ref_dens_liq", "None", ":math:`\text{kg/}\text{m}^3`"
 "Relative density contribution from salt mass fraction", ":math:`\rho_s`", "dens_slope", "None", ":math:`\text{kg/}\text{m}^3`"
 "First density temperature coefficient", ":math:`A`", "dens_param_A", "None", ":math:`\text{K}^{-2}`"
 "Second density temperature coefficient", ":math:`B`", "dens_param_B", "None", ":math:`\text{K}^{-1}`"
 "Third density temperature coefficient", ":math:`C`", "dens_param_C", "None", ":math:`\text{dimensionless}`"
 "First pressure coefficient", ":math:`\alpha`", "ref_pressure_correction", "None", ":math:`\text{dimensionless}`"
 "Second pressure coefficient", ":math:`\beta`", "ref_pressure_slope", "None", ":math:`\text{Pa}^{-1}`"
 "Reference viscosity (at 273 K)", ":math:`\mu_{ref}`", "mu_A", "None", ":math:`\text{kg/}\text{m/}\text{s}`"
 "First viscosity temperature coefficient", ":math:`\mu_B`", "mu_B", "None", ":math:`\text{K}`"
 "Second viscosity temperature coefficient", ":math:`\mu_C`", "mu_C", "None", ":math:`\text{K}`"

**The parameters provided are valid between 0 and 350 °C and up to 600 bar.**


Properties
----------
.. csv-table::
  :header: "Description", "Symbol", "Variable", "Index", "Units"

  "Component mass fraction", ":math:`x_j`", "mass_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`"
  "Mass density of fluid", ":math:`\rho`", "dens_mass_phase", "[p]", ":math:`\text{kg/}\text{m}^3`"
  "Phase volumetric flowrate", ":math:`Q_p`", "flow_vol_phase", "[p]", ":math:`\text{m}^3\text{/s}`"
  "Mass concentration", ":math:`C_j`", "conc_mass_phase_comp", "[p, j]", ":math:`\text{kg/}\text{m}^3`"
  "Dynamic viscosity", ":math:`\mu`", "visc_d_phase", "[p]", ":math:`\text{Pa}\cdotp\text{s}`"
  "Enthalpy flow", ":math:`H`", "enth_flow", "None", ":math:`\text{J/s}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mass fraction", ":math:`x_j = \frac{M_j}{\sum_{j} M_j}`"
   "Mass density", ":math:`\rho = (\rho_{ref} + \rho_s \cdotp \sum_{j} x_j) \cdotp (A \cdotp T^2 + B \cdotp T + C) \cdotp (\alpha + \beta \cdotp P)`"
   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = x_j \cdotp \rho`"
   "Dynamic viscosity", ":math:`\mu = \mu_{ref} \cdotp exp( \frac{\mu_B}{T - \mu_C} )`"
   "Enthalpy flow", ":math:`H = c_p \cdotp \sum_{j} M_j \cdotp (T - 273)`"

Scaling
-------
This coagulation property package includes support for scaling, such as providing
default or calculating scaling factors for almost all variables. The only variables
that do not have scaling factors are the component mass flowrate and the user will
receive a warning if these are not set.

The user can specify the scaling factors for component mass flowrates with the following:

.. testsetup::

  from pyomo.environ import ConcreteModel
  from idaes.core import FlowsheetBlock

.. doctest::

  # relevant imports
  import watertap.property_models.unit_specific.coagulation_prop_pack as props
  from idaes.core.util.scaling import calculate_scaling_factors

  # relevant assignments
  m = ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)
  m.fs.properties = props.CoagulationParameterBlock()

  # set scaling for component mass flowrate
  m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq','H2O'))
  m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq','TDS'))
  m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq','TSS'))
  m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq','Sludge'))

  # calculate scaling factors
  calculate_scaling_factors(m.fs)

The default scaling factors are as follows:

  * 1e-2 for temperature
  * 1e-6 for pressure
  * 1e-3 for mass density
  * 1e3 for dynamic viscosity

Scaling factors for other variables can be calculated based on their relationships
with the user-supplied or default scaling factors.

Reference
---------

Engineering Toolbox. Water - Density, Specific Weight, and
Thermal Expansion Coefficients. (2003)
https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html
[Accessed 02-01-2022]

D.S. Viswananth, G. Natarajan. Data Book on the Viscosity of
Liquids. Hemisphere Publishing Corp. (1989)
