Multi-Component Aqueous Solution (MCAS) Property Package
=========================================================

This property package implements property relationships for an aqueous solution that may contain multiple neutral and/or ionic solutes.

This MCAS property package
   * sets H2O as the solvent;
   * supports multiple solute components including ions and neutral molecules;
   * supports only liquid phase;
   * uses molar flow rate (in mol/s), temperature and pressure as the initial state variables;
   * does not support dynamics.
   
Classes
-------
.. currentmodule:: watertap.property_models.multicomp_aq_sol_prop_pack 

.. autoclass:: MCASParameterBlock
    :members:
    :noindex:

.. autoclass:: MCASParameterData
    :members:
    :noindex:

.. autoclass:: _MCASStateBlock
    :members:
    :noindex:

.. autoclass:: MCASStateBlockData
    :members:
    :noindex:
    
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
   :header: "Description", "Symbol", "Coded Var Name", "Index", "Unit"

   "Component molar flow rate", ":math:`N`", "flow_mol_phase_comp", "[p, j]", ":math:`\text{mol s}^{-1}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"

Calculated Properties
---------------------
.. csv-table::
   :header: "Description", "Symbol", "Coded Var Name", "Index", "Unit", "Calculation Methods"

   "Component mass flow rate", ":math:`M`", "flow_mass_phase_comp", "[p, j]", ":math:`\text{kg s}^{-1}`", ":math:`M=Nm_N`"
   "Component charge-equivalent molar flow rate", ":math:`\tilde{N}`", "flow_equiv_phase_comp", "[p, j]", ":math:`\text{mol s}^{-1}`", ":math:`\tilde{N}=N\left|z\right|`"
   "Component charge-equivalent molar concentration", ":math:`\tilde{n}`", "conc_equiv_phase_comp", "[p, j]", ":math:`\text{mol m}^{-3}`", ":math:`\tilde{n}=n\left|z\right|`"
   "Component mass fraction", ":math:`x`", "mass_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`", ":math:`x_j=\frac{M_j}{\sum_j{M_j}}`"
   "Mass density of aqueous phase", ":math:`\rho`", "dens_mass_phase", "[p]", ":math:`\text{kg m}^{-3}`", ":math:`\rho=1000 \text{kg m}^{-3}` or :math:`\rho=\rho_w + \textbf{f} \left(\sum_{j\in solute}{x_j}, T\right)` :sup:`1`"
   "Phase volumetric flowrate", ":math:`Q`", "flow_vol_phase", "[p]", ":math:`\text{m}^3\text{ } \text{s}^{-1}`", ":math:`Q=\frac{\sum_j{N_j m_{Nj}}}{\rho}`"
   "Total volumetric flowrate", ":math:`Q_{tot}`", "flow_vol", "None", ":math:`\text{m}^3\text{ } \text{s}^{-1}`",":math:`Q_{tot}=\sum_p{Q_p}`" 
   "Component molar concentration", ":math:`n`", "conc_mol_phase_comp", "[p, j]", ":math:`\text{mol m}^{-3}`",":math:`nm_N=m`"
   "Component mass concentration", ":math:`m`", "conc_mass_phase_comp", "[p, j]", ":math:`\text{kg m}^{-3}`",":math:`m=\rho x`"
   "Component molar fraction", ":math:`y`", "mole_frac_phase_comp", "[p, j]", ":math:`\text{dimensionless}`", ":math:`y_j=\frac{N_j}{\sum_j{N_j}}`"
   "Component molality", ":math:`b`", "molality_phase_comp", "[p, j]", ":math:`\text{mol kg}^{-1}`",":math:`b=\frac{N}{N_{H_2O} m_{N\text{H_2O}}}`"
   "Molar volume", ":math:`V_b`", "molar_volume_comp", "['Liq', j]", ":math:`\text{m}^3\text{ } \text{mol}^{-1}`", ":math:`V_b=` input from users"
   "Component diffusivity", ":math:`D`", "diffus_phase_comp", "[p, j]", ":math:`\text{m}^2 \text{ } \text{s}^{-1}`", ":math:`D=` input from users or :math:`D_l\left[ \frac{cm^2}{s} \right] = \frac{13.26\times 10^{-5}}{\left( \mu_w\left[ cP \right] \right)^{1.14}\left( V_b \left[ \frac{cm^3}{mol} \right]\right)^{1.14}}`"
   "Dynamic viscosity", ":math:`\mu`", "visc_d_phase", "[p]", ":math:`\text{Pa s}`", ":math:`\mu=` input from users"
   "Kinematic viscosity", ":math:`\nu`", "visc_k_phase", "[p]", ":math:`\text{m}^2 \text{ s}^{-1}`",":math:`\nu=\mu\rho^{-1}`"
   "Phase osmotic pressure", ":math:`\Pi`", "pressure_osm_phase", "[p]", ":math:`\text{Pa}`",":math:`\Pi=RT\sum_{j\in solute}{n_j}`"
   "Component stokes radius", ":math:`r_h`", "radius_stokes_comp", "[j]", ":math:`\text{m}`", ":math:`r_h=` input from users"
   "Component molecular weight", ":math:`m_N`", "mw_comp", "[j]", ":math:`\text{kg mol}^{-1}`",":math:`m_N=` input from users"
   "Ion component electrical mobility", ":math:`\mu_e`", "elec_mobility_phase_comp", "[p,j]", ":math:`\text{m}^2\text{ }\text{V}^{-1}\text{ }\text{s}^{-1}`", ":math:`\mu_e=` input from users or :math:`\mu_e=\frac{D\left|z\right|F}{RT}`"
   "Ion component transport number", ":math:`t`", "trans_num_phase_comp", "[p, j]", ":math:`\text{dimensionless}`", ":math:`t=` input from users or :math:`t_j=\frac{\left|z_j\right|\mu_{ej} n_j}{\sum_{j\in ion}{\left|z_j\right|\mu_{ej} n_j}}`"
   "Phase equivalent conductivity", ":math:`\Lambda`", "equiv_conductivity_phase", "[p]", ":math:`\text{m}^2 \text{ } \Omega^{-1} \text{ mol}^{-1}`", ":math:`\Lambda=` input from users or :math:`\Lambda=\frac{\sum_{j\in ion}{F\left|z_j\right|\mu_{ej} n_j}}{\sum_{j\in cation}{\left|z_j\right|n_j}}`"
   "Phase electrical conductivity", ":math:`\lambda`", "elec_cond_phase", "[p]", ":math:`\Omega^{-1} \text{ m}^{-1}`", ":math:`\lambda=\Lambda\sum_{j\in cation}{\left|z_j\right|n_j}`"
   "Ion component charge", ":math:`z`", "charge_comp", "[j]", ":math:`\text{dimensionless}`", ":math:`z=` input from users"
   "Component activity coefficient", ":math:`\gamma`", "act_coeff_phase_comp", "[j]", ":math:`\text{dimensionless}`", ":math:`\gamma=` input from users"
   "Dielectric constant", ":math:`\epsilon`", "dielectric_constant", "none", ":math:`\text{dimensionless}`", ":math:`\epsilon=` input from users"
   "Debye-Huckel constant", ":math:`A`", "deby_huckel_constant", "none", ":math:`\text{dimensionless}`", ":math:`A=\frac{\left(2 \pi N_A\right)^{0.5}}{log(10)} \left(\frac{\textbf{e}^2}{4 \pi \epsilon \epsilon_0 kT}\right)^{\frac{3}{2}}`"
   "Ionic Strength", ":math:`I`", "ionic_strength_molal", "none", ":math:`\text{mol kg}^{-1}`", ":math:`I=0.5\sum_{j\in ion}{z_j^2b_j}`"

**Notes** 
   :sup:`1`  :math:`\textbf{f}(\cdot)` refers to empirical correlations of phase or solvent mass density to seawater salinity and temperature following the study of Sharqawy et al. (2010). 

Physical/chemical constants
---------------------------
.. csv-table::
   :header: "Description", "Symbol", "Value", "Unit"
   
   "Idea gas constant", ":math:`R`", "8.3145", ":math:`\text{J mol}^{-1} \text{K}^{-1}`"
   "Faraday constant", ":math:`F`", "96485.33", ":math:`\text{C mol}^{-1}`"
   "Avogadro constant", ":math:`N_A`", "6.022e23", ":math:`\text{dimensionless}`"
   "Boltzmann constant", ":math:`k`", "1.381e-23", ":math:`\text{J K}^{-1}`"
   "Vacuum permittivity", ":math:`\epsilon_0`", "8.854e-12", ":math:`\text{F m}^{-1}`"
   "Elementary charge", ":math:`\textbf{e}`", "1.602e-19", ":math:`\text{C}`"

Scaling
-------
A comprehensive scaling factor calculation method is coded in this property package.  Among the state variables (:math:`N, T, \text{and }  p`), default scaling factors for :math:`T` and :math:`p` were set and do not need users' input, while, for :math:`N`, usually require a user input via an interface. The coding interface to set defalut scaling factor for :math:`N` and call the scaling calculation for other variables is the following. 

.. code-block::

   m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq','{component name}')) 
   # m is the model name, and fs is the instanced flowsheet block of m. 
   calculate_scaling_factors(m)

Users also have the authority to set a scaling factor for non-state variables via the following codes: 

.. code-block::

   import idaes.core.util.scaling as iscale #import the needed utility package
   ...
   iscale.set_scaling_factor(m.fs.properties.{property_name}, 100) 

Proper scaling of variables is, in many cases, crucial to solver's performance in finding an optimal solution of a problem. While designing scaling can have a mathematical sophistication, a general rule is to scale all variables as close to 1 as possible, e.g., in the range of 1e-2 to 1e2. 

   
Reference
---------

M.H. Sharqawy, J.H.L. V, S.M. Zubair, Thermophysical properties of seawater: a review of existing correlations and data, Desalination and Water Treatment. 16 (2010) 354–380. https://doi.org/10.5004/dwt.2010.1079. (2017 corrections provided at http://web.mit.edu/seawater )

Bard, A. J., Faulkner, L. R., & White, H. S. (2022). Electrochemical methods: fundamentals and applications. John Wiley & Sons.

Hayduk, W., & Laudie, H. (1974). Prediction of diffusion coefficients for nonelectrolytes in dilute aqueous solutions. AIChE Journal, 20(3), 611–615. https://doi.org/10.1002/aic.690200329
