Seawater Property Package
=========================

This package implements property relationships for seawater as provided in `Sharqawy et al. (2010) <https://doi.org/10.5004/dwt.2010.1079>`_ and `Nayar et al. (2016) <https://doi.org/10.5004/dwt.2010.1079>`_.

This seawater property package:
   * supports only H2O (solvent) and TDS (solute) components 
   * supports only liquid phase
   * is formulated on a mass basis
   * estimates molar basis properties based on an average molecular weight of sea salt
   * does not support dynamics
   * properties do not incorporate validity ranges for temperature and salinity
   * pressure-dependency of specific enthalpy is incorporated
   * assumes diffusivity of NaCl based on `Bartholomew & Mauter (2019) <https://doi.org/10.1016/j.memsci.2018.11.067>`_

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

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
   "Mass density of seawater", ":math:`\rho`", "dens_mass_phase", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Mass density of pure water", ":math:`\rho_w`", "dens_mass_solvent", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Phase volumetric flowrate", ":math:`Q_p`", "flow_vol_phase", "[p]", ":math:`\text{m}^3\text{/s}`"
   "Volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Mass concentration", ":math:`C_j`", "conc_mass_phase_comp", "[p, j]", ":math:`\text{kg/}\text{m}^3`"
   "Dynamic viscosity", ":math:`\mu`", "visc_d_phase", "[p]", ":math:`\text{Pa}\cdotp\text{s}`"
   "Osmotic coefficient", ":math:`\phi`", "osm_coeff", "None", ":math:`\text{dimensionless}`"
   "Specific enthalpy", ":math:`\widehat{H}`", "enth_mass_phase", "[p]", ":math:`\text{J/kg}`"
   "Enthalpy flow", ":math:`H`", "enth_flow", "None", ":math:`\text{J/s}`"
   "Saturation pressure", ":math:`P_v`", "pressure_sat", "None", ":math:`\text{Pa}`"
   "Specific heat capacity", ":math:`c_p`", "cp_mass_phase", "[p]", ":math:`\text{J/kg/K}`"
   "Thermal conductivity", ":math:`\kappa`", "therm_cond_phase", "[p]", ":math:`\text{W/m/K}`"
   "Latent heat of vaporization", ":math:`h_{vap}`", "dh_vap_mass", "None", ":math:`\text{J/kg}`"
   "Diffusivity", ":math:`D`", "diffus_phase_comp", "[p]", ":math:`\text{m}^2\text{/s}`"
   "Boiling point elevation", ":math:`BPE`", "boiling_point_elevation_phase", "[p]", ":math:`\text{K}`"


**The properties make use of the average molecular weight of sea salt, ≈ 31.40 g/mol, reported in the Reference-Composition Salinity Scale (Millero et al., 2008)  to convert to moles.**

.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mole flowrate", ":math:`N_j`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mole/s}`"
   "Component mole fraction", ":math:`y_j`", "mole_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`" 
   "Molality", ":math:`Cm`", "molality_phase_comp", "['TDS']", ":math:`\text{mole/kg}`"
   "Osmotic pressure", ":math:`\pi`", "pressure_osm_phase", "None", ":math:`\text{Pa}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mass fraction", ":math:`x_j = \frac{M_j}{\sum_{j} M_j}`"
   "Mass density", "Equation 8 in Sharqawy et al. (2010)"
   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = x_j \cdotp \rho`"
   "Dynamic viscosity", "Equations 22 and 23 in Sharqawy et al. (2010)"
   "Osmotic coefficient", "Equation 49 in Sharqawy et al. (2010)"
   "Specific enthalpy", "Equations 25-27 in Nayar et al. (2016)"
   "Enthalpy flow", ":math:`H = \sum_{j} M_j \cdotp \widehat{H}`"
   "Component mole flowrate", ":math:`N_j = \frac{M_j}{MW_j}`"
   "Component mole fraction", ":math:`y_j = \frac{N_j}{\sum_{j} N_j}`"
   "Molality", ":math:`Cm = \frac{x_{TDS}}{(1-x_{TDS}) \cdotp MW_{TDS}}`"
   "Osmotic pressure", ":math:`\pi = \phi \cdotp Cm \cdotp \rho_w \cdotp R \cdotp T` [See note below]"
   "Saturation pressure", "Equations 5 and 6 in Nayar et al. (2016)"
   "Specific heat capacity", "Equation 9 in Sharqawy et al. (2010)"
   "Thermal conductivity", "Equation 13 in Sharqawy et al. (2010)"
   "Latent heat of vaporization", "Equations 37 and 55 in Sharqawy et al. (2010)"
   "Diffusivity", "Equation 6 in Bartholomew et al. (2019)"
   "Boiling point elevation", "Equation 36 in Sharqawy et al. (2010)"



Note: Osmotic pressure calculation (based on equation 48 in Nayar et al. (2016)) uses the density of water as a function of temperature (:math:`\rho_w`) and the ideal gas constant (:math:`R\text{, 8.314 J/mol}\cdotp\text{K}`), in addition to previously defined variables.

Scaling
-------
This seawater property package includes support for scaling, such as providing default or calculating scaling factors for almost all variables. The only variables that do not have scaling factors are the component mass flowrate and the user will receive a warning if these are not set.

The user can specify the scaling factors for component mass flowrates with the following:

.. testsetup::

   from pyomo.environ import ConcreteModel
   from idaes.core import FlowsheetBlock

.. doctest::
   
   # relevant imports
   import watertap.property_models.seawater_prop_pack as props
   from idaes.core.util.scaling import calculate_scaling_factors

   # relevant assignments
   m = ConcreteModel()
   m.fs = FlowsheetBlock(dynamic=False)
   m.fs.properties = props.SeawaterParameterBlock()

   # set scaling for component mass flowrate
   m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq','H2O'))
   m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq','TDS'))

   # calculate scaling factors
   calculate_scaling_factors(m.fs)

The default scaling factors are as follows:

   * 1e-2 for temperature
   * 1e-6 for pressure
   * 1e-3 for mass density
   * 1e3 for dynamic viscosity
   * 1 for the osmotic coefficient
   * 1e-5 for the specific enthalpy
   * 1e-5 for saturation pressure
   * 1e-3 for the specific heat capacity
   * 1 for thermal conductivity
   * 1e-6 for latent heat of vaporization
   * 1e9 for diffusivity
   * 1 for boiling point elevation

Scaling factors for other variables can be calculated based on their relationships with the user-supplied or default scaling factors.
   
References
----------

K.G. Nayar, M.H. Sharqawy, L.D. Banchik, and J.H. Lienhard V, "Thermophysical properties of seawater: A review and new correlations that include pressure dependence,"Desalination, Vol.390, pp.1 - 24, 2016. https://doi.org/10.1016/j.desal.2016.02.024

M.H. Sharqawy, J.H.L. V, S.M. Zubair, Thermophysical properties of seawater: a review of existing correlations and data, Desalination and Water Treatment. 16 (2010) 354–380. https://doi.org/10.5004/dwt.2010.1079. (2017 corrections provided at http://web.mit.edu/seawater)

F.J. Millero, R. Feistel, D.G. Wright, T.J. McDougall, The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Research Part I. 55 (2008) 50–72. https://doi.org/10.1016/j.dsr.2007.10.001.

T.V. Bartholomew, M.S. Mauter, Computational framework for modeling membrane processes without process and solution property simplifications, Journal of Membrane Science. 573 (2019) 682–693. https://doi.org/10.1016/j.memsci.2018.11.067.

