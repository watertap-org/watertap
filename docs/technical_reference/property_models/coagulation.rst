Coagulation Property Package
============================

This package implements property relationships for water density as a function of
temperature and pressure from
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
