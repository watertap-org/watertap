NaCl Property Package with Temperature Dependence
=================================================

This package implements property relationships for an NaCl solution.

This NaCl property package:
   * supports only H2O (solvent) and NaCl (solute) components
   * supports only liquid phase
   * is formulated on a mass basis
   * is intended for non-isothermal applications
   * does not support dynamics


Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Components", ":math:`j`", "['H2O', 'NaCl']"
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
   "Mass density", ":math:`\rho`", "dens_mass_phase", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Phase volumetric flowrate", ":math:`Q_p`", "flow_vol_phase", "[p]", ":math:`\text{m}^3\text{/s}`"
   "Volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Mass concentration", ":math:`C_j`", "conc_mass_phase_comp", "[p, j]", ":math:`\text{kg/}\text{m}^3`"
   "Component mole flowrate", ":math:`N_j`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mole/s}`"
   "Component mole fraction", ":math:`y_j`", "mole_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`" 
   "Molality", ":math:`Cm_{NaCl}`", "molality_comp", "['NaCl']", ":math:`\text{mole/kg}`"
   "Solubility", ":math:`x_{sat}_{NaCl}`", "solubility_comp", "['NaCl']", ":math:`\text{dimensionless}`"
   "Vapor Pressure", ":math:`P_sat`", "pressure_sat", "None", ":math:`\text{Pa}`"
   "Dynamic viscosity", ":math:`\mu`", "visc_d_phase", "[p]", ":math:`\text{Pa/s}`"
   "Thermal Conductivity", ":math:`k`", "therm_cond_phase", "[p]", ":math:`\text{W/m/K}`"
   "Specific Heat Capacity", ":math:`C_p`", "cp_mass_phase", "[p]", ":math:`\text{J/kg/K}`"
   "Solute diffusivity", ":math:`D`", "diffus_phase_comp", "[p,'NaCl']", ":math:`\text{m}^2\text{/s}`"
   "Osmotic coefficient", ":math:`\phi`", "osm_coeff", "None", ":math:`\text{dimensionless}`"
   "Osmotic pressure", ":math:`\pi`", "pressure_osm", "None", ":math:`\text{Pa}`"
   "Specific enthalpy", ":math:`\widehat{H}`", "enth_mass_phase", "[p]", ":math:`\text{J/kg}`"
   "Enthalpy flow", ":math:`H`", "enth_flow", "None", ":math:`\text{J/s}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mass fraction", ":math:`X_j = \frac{M_j}{\sum_{j} M_j}`"
   "Mass density", "Equation 7 in Sparrow (2003)"
   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = X_j \rho`"
   "Component mole flowrate", ":math:`N_j = \frac{M_j}{MW_j}`"
   "Component mole fraction", ":math:`y_j = \frac{N_j}{\sum_{j} N_j}`"
   "Molality", ":math:`Cm_{NaCl} = \frac{x_{NaCl}}{(1-x_{NaCl}) MW_{NaCl}}`"
   "Solubility", "Equation 5 in Sparrow (2003)"
   "Vapor Pressure", "Equation 6 in Sparrow (2003)"
   "Dynamic viscosity", "Regressed from Zaytsev & Aseev (1992)"
   "Thermal Conductivity", "Regressed from Zaytsev & Aseev (1992)"
   "Specific Heat Capacity", "Regressed from Zaytsev & Aseev (1992)"
   "Solute diffusivity", "Regressed from Zaytsev & Aseev (1992)"
   "Osmotic coefficient", "Regressed from Pitzer et. al. (1984)"
   "Osmotic pressure", ":math:`\pi = i \phi Cm_{NaCl} \rho_w R T` [See note 1 below]"
   "Specific enthalpy", "Equation 8 in Sparrow (2003)"
   "Enthalpy flow", ":math:`H = \sum_{j} M_j \widehat{H}`"

Note 1: Osmotic pressure calculation uses the van 't Hoff factor (:math:`i\text{, assumed to be 2}`), density of water (as a function of temperature), gas constant (:math:`R\text{, 8.314 J/mol}\text{K}`) in addition to previously defined variables.


Scaling
-------
This NaCl property package includes support for scaling, such as providing default or calculating scaling factors for almost all variables. The only variables that do not have scaling factors are the component mass flowrate and the user must set them or the user will receive a warning.

The user can specify the scaling factors for component mass flowrates with the following:

The default scaling factors are as follows:

   * 1e-2 for temperature
   * 1e-6 for pressure
   * 1e-3 for mass density
   * 1e3 for dynamic viscosity
   * 1e9 for solute diffusivity
   * 1 for the osmotic coefficient
   * 1e-5 for the specific enthalpy
   * 1e-4 for the specific heat capacity
   * 1e-5 for vapor pressure
   * 1 for thermal conductivity
   * 1 for solubility


The scaling factors for other variables can be calculated based on their relationships with the other variables with the user supplied or default scaling factors.
   
References
----------
.. _Pitzer:

    Pitzer, Kenneth S., J. Christopher Peiper, and R. H. Busey. (1984). Thermodynamic Properties of Aqueous Sodium Chloride Solutions, Journal of Physical and Chemical Reference Data 13, no. 1 , 1–102. `DOI: 10.1063/1.555709 <https://doi.org/10.1063/1.555709>`_

.. _Sharqawy:

   Mostafa H. Sharqawy, John H. Lienhard V & Syed M. Zubair. (2010). Thermophysical properties of seawater: a review of existing correlations and data, Desalination and Water Treatment, 16:1-3, 354-380, `DOI: 10.5004/dwt.2010.1079 <https://doi.org/10.5004/dwt.2010.1079>`_

.. _Sparrow:

    Sparrow, Benjamin S. (2003). Empirical Equations for the Thermodynamic Properties of Aqueous Sodium Chloride, Desalination 159, no. 2, 161–70. `DOI: 10.1016/S0011-9164(03)90068-3 <https://doi.org/10.1016/S0011-9164(03)90068-3>`_

.. _Zaytsev:

    Zaytsev Ivan Dmitrievich & Aseev G. G. (1992). Properties of aqueous solutions of electrolytes, CRC Press.