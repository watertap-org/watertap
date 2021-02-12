Seawater Property Package
=========================

This package implements property relationships for seawater as provided in `Sharqawy et al (2010) <https://doi.org/10.5004/dwt.2010.1079>`_.

This seawater property package:
   * supports only H2O (solvent) and TDS (solute) components 
   * supports only liquid phase
   * is formulated on a mass basis
   * estimates molar basis properties by assuming molecular weight of TDS is equivalent to NaCl
   * does not support dynamics

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indicies"

   "Components", ":math:`j`", "['H2O', 'TDS']"
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
   "Volumetric flowrate", ":math:`Q`", "flow_vol_phase", "[p]", ":math:`\text{m}^3\text{/s}`"
   "Mass concentration", ":math:`C_j`", "conc_mass_phase_comp", "[p, j]", ":math:`\text{kg/}\text{m}^3`"
   "Dynamic viscosity", ":math:`\mu`", "visc_d_phase", "[p]", ":math:`\text{Pa}\cdotp\text{s}`"
   "Osmotic coefficient", ":math:`\phi`", "osm_coeff", "None", ":math:`\text{dimensionless}`"
   "Specific enthalpy", ":math:`\widehat{H}`", "enth_mass_phase", "[p]", ":math:`\text{J/kg}`"
   "Enthalpy flow", ":math:`H`", "enth_flow", "None", ":math:`\text{J/s}`"

**The properties below assume TDS is equivalent to NaCl to convert to moles.**

.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mole flowrate", ":math:`N_j`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mole/s}`"
   "Component mole fraction", ":math:`y_j`", "mole_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`" 
   "Molality", ":math:`Cm_{TDS}`", "molality_comp", "['TDS']", ":math:`\text{mole/kg}`"
   "Osmotic pressure", ":math:`\pi`", "pressure_osm", "None", ":math:`\text{Pa}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mass fraction", ":math:`x_j = \frac{M_j}{\sum_{j} M_j}`"
   "Mass density", "Equation 8 in Sharqawy (2010)"
   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = x_j \cdotp \rho`"
   "Dynamic viscosity", "Equation 22 and 23 in Sharqawy (2010)"
   "Osmotic coefficient", "Equation 49 in Sharqawy (2010)"
   "Specific enthalpy", "Equation 43 and 55 in Sharqawy (2010)"
   "Enthalpy flow", ":math:`H = \sum_{j} M_j \cdotp \widehat{H}`"
   "Component mole flowrate", ":math:`N_j = \frac{M_j}{MW_j}`"
   "Component mole fraction", ":math:`y_j = \frac{N_j}{\sum_{j} N_j}`"
   "Molality", ":math:`Cm_{TDS} = \frac{x_{TDS}}{(1-x_{TDS}) \cdotp MW_{TDS}}`"
   "Osmotic pressure", ":math:`\pi = i \cdotp \phi \cdotp Cm_{TDS} \cdotp \rho_w \cdotp R \cdotp T` [See note below]"

Note: Osmotic pressure calculation uses the van 't Hoff factor (:math:`i\text{, assumed to be 2}`), density of water (:math:`\rho_w\text{, assumed to be 1000 kg/}\text{m}^3`), gas constant (:math:`R\text{, 8.314 J/mol}\cdotp\text{K}`) in addition to previously defined variables.

Scaling
-------
This seawater property package includes support for scaling, such as providing default or calculating scaling factors for almost all variables. The only variables that do not have scaling factors are the component mass flowrate and the user will receive a warning if these are not set.

The user can specify the scaling factors for component mass flowrates with the following:

.. code-block:: python
   
   # relevant imports
   import proteuslib.property_models.seawater_prop_pack as props
   import idaes.core.util.scaling as calculate_scaling_factors

   # relevant assignments
   m = ConcreteModel()
   m.fs = FlowsheetBlock(default={"dynamic": False})
   m.fs.properties = props.SeawaterParameterBlock()

   # set scaling for component mass flowrate
   m.fs.properties.set_default_scaling('flow_mass_comp', 1, index=('Liq','H2O'))
   m.fs.properties.set_default_scaling('flow_mass_comp', 1e2, index=('Liq','TDS'))

   # calculate scaling factors
   calculate_scaling_factors(m.fs)

The default scaling factors are as follows:

   * 1e-2 for temperature
   * 1e-6 for pressure
   * 1e-3 for mass density
   * 1e3 for dynamic viscosity
   * 1 for the osmotic coefficient
   * 1e-5 for the specific enthalpy

The scaling factors for other variables can be calculated based on their relationships with the other variables with the user supplied or default scaling factors.
   
Reference
---------

.. _Sharqawy:

Mostafa H. Sharqawy, John H. Lienhard V & Syed M. Zubair (2010) Thermophysical properties of seawater: a review of existing correlations and data, Desalination and Water Treatment, 16:1-3, 354-380, `DOI: 10.5004/dwt.2010.1079 <https://doi.org/10.5004/dwt.2010.1079>`_
