Multi-Component Aqueous Solution (MCAS) Property Package

============================================

This property package implements property relationships for an aqueous solution that may contain multiple neutral and/or ionic solutes.

This MCAS property package
   * sets H2O as the solvent 
   * supports multiple solute components including ions and neutral molecules.
   * supports only liquid phase
   * uses molar flow rate (in mol/s), temperature and pressure as the initial state variables.  
   * does not support dynamics
   
Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Members"

   "AqueousPhase", ":math:`p`", "{'Liq'}"
   "component_list", ":math:`j`", "{'H2O', solute_list :sup:`1`}"
   "solute_set", ":math:`j`", "{neutral species in solute_list} :sup:`2`"
   "ion_set", ":math:`j`", "{ionic species in solute_list}"
   "cation_set", ":math:`j`", "{cationic species in solute_list}"
   "anion_set", ":math:`j`", "{anionic species in solute_list}"

**Notes** 
   :sup:`1`  solute_list is provided by a necessary configuration to use this property package. 
   :sup:`2` In the implementing codes, the "solute_set" was declared to include only neutral species for current programming convenience; in other places throughout this document, "solute" by itself includes both ionic and neutral species solvated by water.  

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component molar flow rate", ":math:`N`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mols^{-1}}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"

Calculated Properties
---------------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units", "Calculation Methods"

   "Component mass flow rate", ":math:`M`", "flow_mass_phase_comp", "[p, j]", ":math:`\text{kg\ s^{-1}}`", ":math:`M=Nm_N`"
   "Component charge-equivalent molar flow rate", ":math:`\tilde{N}`", "flow_equiv_phase_comp", "[p, j]", ":math:`\text{mol\ s^{-1}}`", ":math:`\tilde{N}=N\left|z\right|`"
   "Component charge-equivalent molar concentration", ":math:`\tilde{n}`", "conc_equiv_phase_comp", "[p, j]", ":math:`\text{mol\ m^{-3}}`", ":math:`\tilde{n}=n\left|z\right|`"
   "Component mass fraction", ":math:`x`", "mass_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`", ":math:`x_j=\frac{M_j}{\sum_j{M_j}}`"
   "Mass density of aqueous phase", ":math:`\rho`", "dens_mass_phase", "[p]", ":math:`\text{kg m^{-3}}`", ":math:`\rho=1000 \text{kg m^{-3}}` or :math:`\rho=\rho_w \textbf{f} \left(\sum_{j\in solute}{x_j}, T\right)`"
   "Mass density of solvent water", ":math:`\rho_w`", "dens_mass_w_phase", "[p]", ":math:`\text{kg m^{-3}}`",":math:`\rho_w=\textbf{f}\left(T\right)`"
   "Phase volumetric flowrate", ":math:`Q`", "flow_vol_phase", "[p]", ":math:`\text{m^3 s^{-1}}`", ":math:`Q=\frac{\sum_j{N_j m_{Nj}}}{\rho}`"
   "Total volumetric flowrate", ":math:`Q_tot`", "flow_vol", "None", ":math:`\text{m^3 s^{-1}}`",":math:`Q_tot=\sum_p{Q_p}`" 
   "Component molar concentration", ":math:`n`", "conc_mol_phase_comp", "[p, j]", ":math:`\text{mol m^{-3}}`",":math:`nm_N=m`"
   "Component mass concentration", ":math:`m`", "conc_mass_phase_comp", "[p, j]", ":math:`\text{kg m^{-3}}`",":math:`m=\rho x`"
   "Component molar fraction", ":math:`y`", "mole_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`", ":math:`y_j=\frac{N_j}{\sum_j{N_j}}`"
   "Component molality", ":math:`b`", "molality_phase_comp", "[p, j]", ":math:`\text{mol kg^{-1}}`",":math:`b=N\left(N_{H_2O}m_N _{\text{H_2O}}\right)^{-1}`"
   "Component diffusivity", ":math:`D`", "diffus_phase_comp", "[p, j]", ":math:`\text{m^2\ s^{-1}}`", ":math:`D=` input from users"
   "Dynamic viscosity", ":math:`\mu`", "visc_d_phase", "[p]", ":math:`\text{Pa s}`", ":math:`\mu=` input from users"
   "Kinematic viscosity", ":math:`\nu`", "visc_k_phase", "[p]", ":math:`\text{m^2\ s^{-1}}`",":math:`\nu=\mu\rho^{-1}`"
   "Phase osmotic pressure", ":math:`\Pi`", "pressure_osm_phase", "None", ":math:`\text{Pa}`",":math:`\Pi=RT\sum_{j\in solute}{n_j}`"
   "Component stokes radius", ":math:`r_h`", "radius_stokes_comp", "[j]", ":math:`\text{m}`", ":math:`r_h=` input from users"
   "Component molecular weight", ":math:`m_N`", "mw_comp", "[j]", ":math:`\text{kg mol^{-1}}`",":math:`m_N=` input from users"
   "Ion component electrical mobility", ":math:`\mu_e`", "elec_mobility_phase_comp", "[p,j]", ":math:`\text{m^2\ V^{-1}\ s^{-1}}`", ":math:`\mu_e=` input from users or :math:`\mu_e=\frac{D\left|z\right|F}{RT}`"
   "Ion component transport number", ":math:`t`", "trans_num_phase_comp", "[p, j]", ":math:`\text{dimensionless}`", ":math:`t=` input from users or :math:`t_j=\frac{\left|z_j\right|\mu_e _j n_j}{\sum_{j\in ion}{\left|z_j\right|\mu_e _j n_j}}`"
   "Phase equivalent conductivity", ":math:`\Lambda`", "equiv_conductivity_phase", "[p]", ":math:`\text{m^2\ \Omega ^{-1}\  mol^{-1}}`", ":math:`\Lambda=` input from users or :math:`\Lambda=\frac{\sum_{j\in ion}{F\left|z_j\right|\mu_e_j n_j}}{\sum_{j\in cation}{\left|z_j\right|n_j}}`"
   "Phase electrical conductivity", ":math:`\lambda`", "elec_cond_phase", "[p]", ":math:`\text{\Omega^{-1}\ m^{-1}}`", ":math:`\lambda=\Lambda\sum_{j\in cation}{\left|z_j\right|n_j}`"
   "Ion component charge", ":math:`z`", "charge_comp", "[j]", ":math:`\text{dimensionless}`", ":math:`z=` input from users"
   "Component activity coefficient", ":math:`\gamma`", "act_coeff_phase_comp", "[j]", ":math:`\text{dimensionless}`", ":math:`\gamma=` input from users"
   "Dielectric constant", ":math:`\epsilon`", "dielectric_constant", "none", ":math:`\text{dimensionless}`", ":math:`\epsilon=` input from users"
   "Debye-Huckel constant", ":math:`A`", "deby_huckel_constant", "none", ":math:`\text{dimensionless}`", ":math:`A=\frac{\left(2 \pi N_A\right)^{0.5}}{log(10)} \left(\frac{\textbf{e}^2}{4 \pi \epsilon \epsilon_0 kT}\right)^{\left(\frac{3}{2}\right)}`"
   "Ionic Strength", ":math:`I`", "ionic_strength_molal", "none", ":math:`\text{mol kg^{-1}}`", ":math:`I=0.5\sum_{j\in solute}{z_j^2b_j}`"




   "Component mass fraction", ":math:`x = \frac{M}{\sum_{j} M}`"
   "Mass density", "Equation 8 in Sharqawy et al. (2010)"
   "Volumetric flow rate", ":math:`Q = \frac{\sum_{j} M}{\rho}`"
   "Mass concentration", ":math:`C = x \cdotp \rho`"
   "Dynamic viscosity", "Equations 22 and 23 in Sharqawy et al. (2010)"
   "Osmotic coefficient", "Equation 49 in Sharqawy et al. (2010)"
   "Specific enthalpy", "Equations 43 and 55 in Sharqawy et al. (2010)"
   "Enthalpy flow", ":math:`H = \sum_{j} M \cdotp \widehat{H}`"
   "Component mole flow rate", ":math:`N = \frac{M}{MW}`"
   "Component mole fraction", ":math:`y = \frac{N}{\sum_{j} N}`"
   "Molality", ":math:`Cm = \frac{x_{TDS}}{(1-x_{TDS}) \cdotp MW_{TDS}}`"
   "Osmotic pressure", ":math:`\pi = \phi \cdotp Cm \cdotp \rho_w \cdotp R \cdotp T` [See note below]"
   "Saturation pressure", "Equations 5 and 6 in Nayar et al. (2016)"
   "Specific heat capacity", "Equation 9 in Sharqawy et al. (2010)"
   "Thermal conductivity", "Equation 13 in Sharqawy et al. (2010)"
   "Latent heat of vaporization", "Equations 37 and 55 in Sharqawy et al. (2010)"
   "Diffusivity", "Equation 6 in Bartholomew et al. (2019)"



Note: Osmotic pressure calculation (based on equation 48 in Nayar et al. (2016)) uses the density of water as a function of temperature (:math:`\rho_w`) and the ideal gas constant (:math:`R\text{, 8.314 J/mol}\cdotp\text{K}`), in addition to previously defined variables.

Scaling
-------
This seawater property package includes support for scaling, such as providing default or calculating scaling factors for almost all variables. The only variables that do not have scaling factors are the component mass flow rate and the user will receive a warning if these are not set.

The user can specify the scaling factors for component mass flow rates with the following:

.. testsetup::

   from pyomo.environ import ConcreteModel
   from idaes.core import FlowsheetBlock

.. doctest::
   
   # relevant imports
   import watertap.property_models.seawater_prop_pack as props
   from idaes.core.util.scaling import calculate_scaling_factors

   # relevant assignments
   m = ConcreteModel()
   m.fs = FlowsheetBlock(default={"dynamic": False})
   m.fs.properties = props.SeawaterParameterBlock()

   # set scaling for component mass flow rate
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

Scaling factors for other variables can be calculated based on their relationships with the user-supplied or default scaling factors.
   
Reference
---------

K.G.Nayar, M.H.Sharqawy, L.D.Banchik, and J.H.Lienhard V, "Thermophysical properties of seawater: A review and new correlations that include pressure dependence,"Desalination, Vol.390, pp.1 - 24, 2016. https://doi.org/10.1016/j.desal.2016.02.024

M.H. Sharqawy, J.H.L. V, S.M. Zubair, Thermophysical properties of seawater: a review of existing correlations and data, Desalination and Water Treatment. 16 (2010) 354–380. https://doi.org/10.5004/dwt.2010.1079. (2017 corrections provided at http://web.mit.edu/seawater )

F.J. Millero, R. Feistel, D.G. Wright, T.J. McDougall, The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Research Part I. 55 (2008) 50–72. https://doi.org/10.1016/j.dsr.2007.10.001.

T.V. Bartholomew, M.S. Mauter, Computational framework for modeling membrane processes without process and solution property simplifications, Journal of Membrane Science. 573 (2019) 682–693. https://doi.org/10.1016/j.memsci.2018.11.067.

