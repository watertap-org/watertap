Water Property Package
======================

This package implements property relationships for pure water.

This water property package:
   * supports only H2O
   * supports liquid and vapor phases
   * is formulated on a mass basis
   * does not support dynamics
   * pressure-dependency of specific enthalpy is incorporated

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Components", ":math:`j`", "['H2O']"
   "Phases", ":math:`p`", "['Liq', 'Vap']"

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

   "Mass density of pure water", ":math:`\rho`", "dens_mass_phase", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Phase volumetric flowrate", ":math:`Q_p`", "flow_vol_phase", "[p]", ":math:`\text{m}^3\text{/s}`"
   "Volumetric flowrate", ":math:`Q`", "flow_vol", "None", ":math:`\text{m}^3\text{/s}`"
   "Specific enthalpy", ":math:`\widehat{H}`", "enth_mass_phase", "[p]", ":math:`\text{J/kg}`"
   "Enthalpy flow", ":math:`H`", "enth_flow_phase", "[p]", ":math:`\text{J/s}`"
   "Saturation pressure", ":math:`P_v`", "pressure_sat", "None", ":math:`\text{Pa}`"
   "Specific heat capacity", ":math:`c_p`", "cp_mass_phase", "[p]", ":math:`\text{J/kg/K}`"
   "Latent heat of vaporization", ":math:`h_{vap}`", "dh_vap_mass", "None", ":math:`\text{J/kg}`"
   "Component mole flowrate", ":math:`N_j`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mole/s}`"
   "Component mole fraction", ":math:`y_j`", "mole_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`"
   "Dynamic viscosity", ":math:`\mu`", "visc_d_phase", "[p]", ":math:`\text{Pa}\cdotp\text{s}`" 
   "Thermal conductivity", ":math:`\kappa`", "therm_cond_phase", "[p]", ":math:`\text{W/m/K}`"
   
Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mass fraction", ":math:`x_j = \frac{M_j}{\sum_{j} M_j}`"
   "Mass density of liquid", "Equation 8 in Sharqawy et al. (2010)"
   "Mass density of vapor*", ":math:`\rho = \frac{Pm}{nRT}`"
   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = x_j \cdotp \rho`"
   "Specific enthalpy of liquid", "Equations 25-27 in Nayar et al. (2016)"
   "Specific enthalpy of vapor", "Equations 25-27 in Nayar et al. (2016) :math:`+ h_{vap}`"
   "Enthalpy flow", ":math:`H = \sum_{j} M_j \cdotp \widehat{H}`"
   "Component mole flowrate", ":math:`N_j = \frac{M_j}{MW_j}`"
   "Component mole fraction", ":math:`y_j = \frac{N_j}{\sum_{j} N_j}`"
   "Saturation pressure", "Equation 6 in Nayar et al. (2016)"
   "Specific heat capacity of liquid", "Equation 9 in Sharqawy et al. (2010)"
   "Specific heat capacity of vapor", "Shomate equation from NIST WebBook"
   "Latent heat of vaporization", "Equations 37 and 55 in Sharqawy et al. (2010)"
   "Dynamic viscosity", "Equations 22 and 23 in Sharqawy et al. (2010)"
   "Thermal conductivity", "Equation 13 in Sharqawy et al. (2010)"

\* Derived from the ideal gas law


Scaling
-------
This water property package includes support for scaling, such as providing default or calculating scaling factors for almost all variables. The only variables that do not have scaling factors are the component mass flowrate and the user will receive a warning if these are not set.

The user can specify the scaling factors for component mass flowrates with the following:

.. testsetup::

   from pyomo.environ import ConcreteModel
   from idaes.core import FlowsheetBlock

.. doctest::
   
   # relevant imports
   import watertap.property_models.water_prop_pack as props
   from idaes.core.util.scaling import calculate_scaling_factors

   # relevant assignments
   m = ConcreteModel()
   m.fs = FlowsheetBlock(dynamic=False)
   m.fs.properties = props.WaterParameterBlock()

   # set scaling for component mass flowrate
   m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq','H2O'))

   # calculate scaling factors
   calculate_scaling_factors(m.fs)

The default scaling factors are as follows:

   * 1e-2 for temperature
   * 1e-5 for pressure
   * 1e-3 for liquid mass density
   * 1 for vapor mass density
   * 1e-5 for the liquid specific enthalpy
   * 1e-6 for the vapor specific enthalpy
   * 1e-5 for saturation pressure
   * 1e-3 for the liquid specific heat capacity
   * 1e-3 for the vapor specific heat capacity
   * 1e-6 for latent heat of vaporization
   * 1e3 for the dynamic viscosity
   * 1 for the thermal conductivity
  
Scaling factors for other variables can be calculated based on their relationships with the user-supplied or default scaling factors.
   
References
----------

K.G. Nayar, M.H. Sharqawy, L.D. Banchik, and J.H. Lienhard V, "Thermophysical properties of seawater: A review and new correlations that include pressure dependence,"Desalination, Vol.390, pp.1 - 24, 2016. https://doi.org/10.1016/j.desal.2016.02.024

M.H. Sharqawy, J.H.L. V, S.M. Zubair, Thermophysical properties of seawater: a review of existing correlations and data, Desalination and Water Treatment. 16 (2010) 354–380. https://doi.org/10.5004/dwt.2010.1079. (2017 corrections provided at http://web.mit.edu/seawater)

F.J. Millero, R. Feistel, D.G. Wright, T.J. McDougall, The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Research Part I. 55 (2008) 50–72. https://doi.org/10.1016/j.dsr.2007.10.001.

T.V. Bartholomew, M.S. Mauter, Computational framework for modeling membrane processes without process and solution property simplifications, Journal of Membrane Science. 573 (2019) 682–693. https://doi.org/10.1016/j.memsci.2018.11.067.

Water Gas Phase Thermochemistry Data, National Institute of Standards and Technology, 2021, https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&amp;Mask=1#Thermo-Gas.

