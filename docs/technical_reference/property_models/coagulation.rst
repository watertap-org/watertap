Coagulation Property Package
============================

This package implements property relationships for water density as a function of
temperature, pressure, and mass fraction of suspended/dissolved solids from
`Engineering Toolbox. (2003) <https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html>`_
and water viscosity as a function of temperature from
`D.S. Viswananth, G. Natarajan. (1989) <https://www.osti.gov/biblio/6562161>`_.

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
   "Mass density", "blah"
   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = x_j \cdotp \rho`"
   "Dynamic viscosity", "blah"
   "Enthalpy flow", ":math:`H = \sum_{j} M_j \cdotp \widehat{H}`"
