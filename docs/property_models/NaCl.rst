NaCl Property Package
=====================

This package implements property relationships for a NaCl solution as provided in `Bartholomew and Mauter (2019) <https://doi.org/10.1016/j.memsci.2018.11.067>`_.

The NaCl property package has only H2O (solvent) and NaCl (solute) components and only liquid phase. The package is based on mass flowrate and is in SI units. The package is currently not compatible for dynamics.

For this package, the NaCl solution properties are a function of salinity only. The package is based on regressions of experimental data at a temperature of 25 C. Using this package for higher or lower temperatures will result in larger errors. 

This package includes one temperature dependent property, specific enthalpy. The specific enthalpy is based on relationships developed for seawater instead of a NaCl solution, but the deviation is expected to be negiligible for isothermal applications (the recommended use for this property package).

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mass flowrate", ":math:`M_j`", "flow_mass_comp[j]", "j in [H2O, TDS]", ":math:`\text{kg/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"

Properties
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mass fraction", ":math:`X_j`", "mass_frac_comp[j]", "j in [H2O, TDS]", ":math:`\text{dimensionless}`"
   "Mass density", ":math:`\rho`", "dens_mass", "None", ":math:`\text{kg/}\text{m}^3`"
   "Volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Mass concentration", ":math:`C_j`", "conc_mass_comp", "j in [H2O, TDS]", ":math:`\text{kg/}\text{m}^3`"
   "Dynamic viscosity", ":math:`\mu`", "visc_d", "None", ":math:`\text{Pa}\cdotp\text{s}`"
   "Solute diffusivity", ":math:`D`", "diffus", "None", ":math:`\text{m}^2\text{/s}`"
   "Osmotic coefficient", ":math:`\phi`", "osm_coeff", "None", ":math:`\text{dimensionless}`"
   "Osmotic pressure", ":math:`\pi`", "pressure_osm", "None", ":math:`\text{Pa}`"
   "Specific enthalpy", ":math:`\widehat{H}`", "enth_mass", "None", ":math:`\text{J/kg}`"
   "Enthalpy flow", ":math:`H`", "enth_flow", "None", ":math:`\text{J/s}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mass fraction", ":math:`X_j = \frac{M_j}{\sum_{j} M_j}`"
   "Mass density", "Equation 4 in Bartholomew & Mauter (2019)"
   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = X_j \cdotp \rho`"
   "Dynamic viscosity", "Equation 5 in Bartholomew & Mauter (2019)"
   "Solute diffusivity", "Equation 6 in Bartholomew & Mauter (2019)"
   "Osmotic coefficient", "Equation 3b in Bartholomew & Mauter (2019)"
   "Osmotic pressure", ":math:`\pi = i \cdotp \phi \cdotp C_m \cdotp \rho_w \cdotp R \cdotp T` [See note 1 below]"
   "Specific enthalpy", "Equation 43 and 55 in Sharqawy (2010) [See note 2 below]"
   "Enthalpy flow", ":math:`H = \sum_{j} M_j \cdotp \widehat{H}`"

Note 1: Osmotic pressure calculation uses the Van't Hoff factor (:math:`i`), molality (:math:`C_m`), density of water (:math:`\rho_w`), gas constant (:math:`R`) in addition to previously defined variables.

Note 2: Specific enthalpy calculation is based on seawater relationships in Sharqawy (2010).

Scaling
-------
This seawater property package includes support for scaling, such as providing default or calculating scaling factors for almost all variables. The only variables that do not have scaling factors are the component mass flowrate and the user must set them or the user will receive a warning.

The user can specify the scaling factors for component mass flowrates with the following:

.. code-block:: python
   
   # relevant imports
   import proteuslib.property_models.NaCl_prop_pack as props
   import idaes.core.util.scaling as calculate_scaling_factors

   # relevant assignments
   m = ConcreteModel()
   m.fs = FlowsheetBlock(default={"dynamic": False})
   m.fs.properties = props.SeawaterParameterBlock()

   # set scaling for component mass flowrate
   m.fs.properties.set_default_scaling('flow_mass_comp', 1, index='H2O')
   m.fs.properties.set_default_scaling('flow_mass_comp', 1e2, index='TDS')

   # calculate scaling factors
   calculate_scaling_factors(m.fs)

The default scaling factors are as follows:

   * 1e-2 for temperature
   * 1e-6 for pressure
   * 1e-3 for mass density
   * 1e3 for dynamic viscosity
   * 1e9 for solute diffusivity
   * 1 for the osmotic coefficient
   * 1e-5 for the specific enthalpy

The scaling factors for other variables can be calculated based on their relationships with the other variables with the user supplied or default scaling factors.
   
References
----------

.. _Bartholomew:

   Timothy V. Bartholomew, Meagan S. Mauter (2019) Computational framework for modeling membrane processes without process and solution property simplifications, Journal of Membrane Science, 573, 682-693, `DOI: 10.1016/j.memsci.2018.11.067 <https://doi.org/10.1016/j.memsci.2018.11.067>`_

.. _Sharqawy:

   Mostafa H. Sharqawy, John H. Lienhard V & Syed M. Zubair (2010) Thermophysical properties of seawater: a review of existing correlations and data, Desalination and Water Treatment, 16:1-3, 354-380, `DOI: 10.5004/dwt.2010.1079 <https://doi.org/10.5004/dwt.2010.1079>`_
