Granular Activated Carbon (GAC)
===============================
.. note::

    The GAC model is an ongoing project and subject to change, as is this documentation.

This is a empirical-performance-based granular activated carbon (GAC) model that works under the following criteria and assumptions
   * supports a single liquid phase only
   * supports a single solute and single solvent (water) only
   * supports steady-state only
   * assumes isothermal conditions
   * model performance is independent of a gravity-fed or pressurized GAC unit, therefore assumes isobaric

.. index::
   pair: watertap.unit_models.gac;gac

.. currentmodule:: watertap.unit_models.gac

Introduction
------------
The implemented model for estimating GAC performance is adapted from a simplified model originally presented in Hand, 1984 and
further elaborated in Crittenden, 2012. This formulation is denoted as the constant pattern homogeneous surface diffusion model (CPHSDM).
As a GAC system is operated as a batch process, a mass transfer zone (MTZ) is formed in the bed where a concentration profile, or
breakthrough curve, developeds as a function of the adsorption properties. This MTZ is bounded by saturated GAC upstream and
fresh GAC downstream. The CPHSDM is valid under the assumption that the shape of the MTZ is constant as it travels through the
bed and a constant pattern solution (CPS) may be determined. The CPS is calculated through a multistep procedure utilizing common
dimensionless groups applied in polynomial fits to determine performance. Therefore, coefficients used in the polynomial must be
derived from experimental dataof the intended system to produce valid results. Coefficients for common compounds treated by GAC may be
found in both Hand, 1984 and Crittenden, 2012. The model is estimated to have within 10% error and therefore may be applied to bed lengths
shorter than the minimum length determined by the CPHSDM within the error threshold (or lengths greater than the minimum).

The batch operation results of the CPS are converted to approximate steady state results for intuitive use of the model for
flowsheet purposes. A description of the transformation is provided in Figure 1.

.. figure:: ../../_static/unit_models/gac.png
    :width: 1200
    :align: center

    Figure 1. (a) The mass transfer zone movement as a constant pattern for an arbitrary time of operation where the entire MTZ
    is contained within the bed length. The dimensionless concentration of the adsorbent in both the liquid and adsorbate phase is
    variable in the range of the MTZ as plotted with respect to the dimensionless length of the bed. (b) The breakthrough curve is
    shown as a function of the elapsed time in operation of the bed. The dashed line indicates a breakthrough time where the bed operation will
    be stopped (and the partially saturated GAC will be replaced) when the effluent concentration is 50% of the inlet. (c) For the breakthrough
    time in b, some of the MTZ is still contained within the bed and therefore the GAC particles in this zone are not fully saturated. The
    trapezoid rule for integration is used to approximate the degree of saturation of the entire bed. This degree of saturation corresponds
    to a specific mass of contaminant adsorbed during the operation time and is used to transform the process to a steady state approximation.

Degrees of Freedom
------------------
In the default configuration of the GAC unit model there are 18 degrees of freedom in addition to the inlet state variables
(i.e. temperature, pressure, component flowrates) that should be fixed for the model to be fully specified.
In association with using the Freundlich adsorption isotherm and empirical model, the following 9 variables are almost always fixed
and may be derived from experimental data:

   * Freundlich isotherm k parameter
   * Freundlich isotherm 1/n parameter
   * Stanton number equation parameters :math:`\text{a}_0` and :math:`\text{a}_1`
   * throughput ratio equation parameters :math:`\text{b}_0`, :math:`\text{b}_1`, :math:`\text{b}_2`, :math:`\text{b}_3` and :math:`\text{b}_4`

Additionally, the following 9 variables are traditionally fixed:

   * the target dimensionless concentration *or* the dimensionless concentration at the time of replacement
   * empty bed contact time
   * bed voidage
   * superficial velocity
   * GAC particle porosity
   * GAC particle apparent density *or* GAC particle solid density
   * GAC particle diameter
   * liquid phase film transfer coefficient
   * surface diffusion coefficient

When setting the configuration options to calculate the liquid phase film transfer coefficient and surface diffusion coefficient,
these respective variables are no longer specified and 4 newly introduced variables must be fixed. This is a net result of 20 degrees of freedom.
Newly utilized variables that must be fixed include:

   * molal volume of the solute
   * GAC particle sphericity
   * tortuosity of the path that the adsorbate must take as compared to the radius
   * surface-to-pore diffusion flux ratio

Model Structure
------------------
The GAC model consists of 1 ControlVolume0DBlock (process_flow) for the process flow of the water treatment train.
The process flow includes 2 StateBlocks (inlet and outlet) which are used for mass and momentum balances.
It also includes 1 StateBlock (adsorbed) for the solute that is adsorbed into the GAC parties.
The material removed in the adsorbed state block is simulated as liquid phase solute but should be interpreted as solute that has adsorbed
into the solid phase GAC particles. The steady state mass removal and replacement rate of the GAC itself is provided as a unit model
variable and excluded from flowsheet material balances. Therefore, GAC is never included as a component in property package specifications.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'Solute']*"

\*Solute is provided in the imported property model.

.. _GAC_variables:

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Freundlich isotherm k parameter", ":math:`k`", "freund_k", "None", ":math:`\left(\text{m}^3\text{/kg}\right)^\left( \frac{1}{n} \right)`"
   "Freundlich isotherm 1/n parameter", ":math:`\frac{1}{n}`", "freund_ninv", "None", ":math:`\text{dimensionless}`"
   "Equilibrium concentration of adsorbed phase with liquid phase", ":math:`q_e`", "equil_conc", "None", ":math:`\left( \text{kg}_\text{adsorbate}\text{/kg}_\text{adsorbent} \right)`"
   "Mass of contaminant absorbed at the time of replacement", ":math:`M_{solute}`", "mass_adsorbed", "None", ":math:`\text{kg}`"
   "Mass of contaminant adsorbed if fully saturated", ":math:`M_{solute\text{,}e}`", "mass_adsorbed_saturated", "None", ":math:`\text{kg}`"
   "Adsorber bed void fraction", ":math:`\epsilon`", "bed_voidage", "None", ":math:`\text{dimensionless}`"
   "Adsorber bed volume", ":math:`V`", "bed_volume", "None", ":math:`\text{m}^3`"
   "Adsorber bed area", ":math:`A`", "bed_area", "None", ":math:`\text{m}^2`"
   "Adsorber bed length", ":math:`L`", "bed_length", "None", ":math:`\text{m}`"
   "Mass of fresh GAC in the bed", ":math:`M`", "bed_mass_gac", "None", ":math:`\text{kg}`"
   "Superficial velocity", ":math:`v_s`", "velocity_sup", "None", ":math:`\text{m/s}`"
   "Interstitial velocity", ":math:`v_i`", "velocity_int", "None", ":math:`\text{m/s}`"
   "GAC particle porosity", ":math:`\epsilon_p`", "particle_porosity", "None", ":math:`\text{dimensionless}`"
   "GAC particle apparent density", ":math:`\rho_a`", "particle_dens_app", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle bulk density", ":math:`\rho_b`", "particle_dens_bulk", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle solid density", ":math:`\rho_s`", "particle_dens_sol", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle diameter", ":math:`d_p`", "particle_dia", "None", ":math:`\text{m}`"
   "Average dimensionless concentration in the duration of operation", ":math:`\bar{\left( \frac{C}{C_{0}} \right)}`", "conc_ratio_avg", "None", ":math:`\text{dimensionless}`"
   "Dimensionless concentration of the effluent at the time of replacement", ":math:`\left( \frac{C}{C_{0}} \right)|_{z=L}`", "conc_ratio_avg", "None", ":math:`\text{dimensionless}`"
   "Approximate saturation of the GAC in the bed at the time of replacement", ":math:`\bar{\left( \frac{q}{q_{e}} \right)}`", "gac_saturation_replace", "None", ":math:`\text{dimensionless}`"
   "Empty bed contact time", ":math:`EBCT`", "ebct", "None", ":math:`\text{s}`"
   "Mass throughput ratio", ":math:`T`", "mass_throughput", "None", ":math:`\text{dimensionless}`"
   "Residence time", ":math:`\tau`", "res_time", "None", ":math:`\text{s}`"
   "Elapsed operation time", ":math:`t`", "elap_time", "None", ":math:`\text{s}`"
   "Steady state GAC replacement rate", ":math:`\dot{m}_{GAC}`", "gac_mass_replace_rate", "None", ":math:`\text{s}`"
   "Liquid phase film transfer coefficient", ":math:`k_f`", "kf", "None", ":math:`\text{m/s}`"
   "Surface diffusion coefficient", ":math:`D_s`", "ds", "None", ":math:`\text{m}^2\text{/s}`"
   "Solute distribution parameter", ":math:`D_g`", "dg", "None", ":math:`\text{dimensionless}`"
   "Biot number", ":math:`Bi`", "N_Bi", "None", ":math:`\text{dimensionless}`"
   "Minimum Stanton number for CPS", ":math:`St_{min}`", "min_N_St", "None", ":math:`\text{dimensionless}`"
   "Minimum empty bed contact time for CPS", ":math:`EBCT_{min}`", "min_ebct", "None", ":math:`\text{s}`"
   "Minimum residence time for CPS", ":math:`\tau_{min}`", "min_res_time", "None", ":math:`\text{s}`"
   "Minimum elapsed operation time for CPS", ":math:`t_{min}`", "min_elap_time", "None", ":math:`\text{s}`"
   "Mass throughput ratio for the upstream edge of the MTZ", ":math:`T|\left( \frac{C}{C_{0}} = 0.95 \right)`", "mass_throughput_mtz_upstream", "None", ":math:`\text{dimensionless}`"
   "Empty bed contact time of the partial MTZ at the time of replacement", ":math:`EBCT_{MTZ}`", "ebct_mtz_replace", "None", ":math:`\text{s}`"
   "Length of the partial MTZ at the time of replacement", ":math:`L_{MTZ}`", "length_mtz_replace", "None", ":math:`\text{m}`"
   "Stanton equation parameter 0", ":math:`a_0`", "a0", "None", ":math:`\text{dimensionless}`"
   "Stanton equation parameter 1", ":math:`a_1`", "a1", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 0", ":math:`b_0`", "b0", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 1", ":math:`b_1`", "b1", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 2", ":math:`b_2`", "b2", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 3", ":math:`b_3`", "b3", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 4", ":math:`b_4`", "b4", "None", ":math:`\text{dimensionless}`"

The following variables are only built when specific configuration options are selected.

if ``film_transfer_coefficient_type`` or ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Molecular diffusion coefficient", ":math:`D_l`", "diffus_liq", "None", ":math:`\text{m}^2\text{/s}`"
   "Molal volume", ":math:`V_b`", "molal_volume", "None", ":math:`\text{m}^3\text{/mol}`"

if ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Reynolds number", ":math:`Re`", "N_Re", "None", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc`", "N_Sc", "None", ":math:`\text{dimensionless}`"
   "GAC particle sphericity", ":math:`\phi`", "sphericity", "None", ":math:`\text{dimensionless}`"

if ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Tortuosity of the path that the adsorbate must take as compared to the radius", ":math:`\tau_p`", "tort", "None", ":math:`\text{dimensionless}`"
   "Surface-to-pore diffusion flux ratio", ":math:`SPDFR`", "spdfr", "None", ":math:`\text{dimensionless}`"

.. _GAC_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "Equilibrium concentration", ":math:`q_e = kC_0^{1/n}`"
   "Solute distribution parameter", ":math:`D_g=\frac{\rho_aq_e\left( 1-\epsilon \right)}{\epsilon C_0}`"
   "Biot number", ":math:`Bi=\frac{k_fd_p\left( 1-\epsilon \right)}{2D_sD_g\epsilon}`"
   "Minimum Stanton number to achieve CPS", ":math:`St_{min}=a_0Bi+a_1`"
   "Minimum EBCT to achieve CPS", ":math:`EBCT_{min}=\frac{St_{min}d_p}{2k_f\left( 1-\epsilon \right)}`"
   "Minimum residence time to achieve CPS", ":math:`\tau_{min}=EBCT_{min}\epsilon`"
   "Residence time", ":math:`\tau=EBCT\epsilon`"
   "Mass throughput ratio", ":math:`T=b_0+b_1\left( \frac{C}{C_0} \right)^{b_2}+\frac{b_3}{1.01-\left( \frac{C}{C_0} \right)^{b_4}}`"
   "Minimum elapsed operation time to achieve CPS", ":math:`t_{min}=\tau_{min}\left( D_g+1 \right)T`"
   "Elapsed operation time", ":math:`t=t_{min}+\left( \tau-\tau_{min} \right)\left( D_g+1 \right)`"
   "Density relation to bed voidage", ":math:`\epsilon=1-\frac{\rho_b}{\rho_a}`"
   "Density relation to particle porosity", ":math:`\epsilon_p=1-\frac{\rho_a}{\rho_s}`"
   "Steady state GAC replacement rate", ":math:`\dot{m}_{GAC}=\frac{M_{GAC}}{t}`"
   "Adsorber bed volume", ":math:`EBCT=\frac{V}{Q}`"
   "Mass of GAC in a fresh adsorber bed", ":math:`M_{GAC}=V\rho_b`"
   "Velocity relationship", ":math:`v_i=\frac{v_s}{\epsilon}`"
   "Adsorbed bed area", ":math:`A=\frac{Q}{v_s}`"
   "Adsorbed bed length", ":math:`EBCT=\frac{L}{v_s}`"
   "Mass of solute adsorbed if the bed was fully saturated", ":math:`M_{solute\text{,}e}=q_eM_{GAC}`"
   "Saturation fraction of the bed at the time of replacement", ":math:`\bar{\left( \frac{q}{q_{e}} \right)}=\frac{M_{solute}}{M_{solute\text{,}e}}`"
   "Mass throughput ratio of the upstream edge of the MTZ", ":math:`T_{MTZ}=b_0+b_1\left( \frac{C}{C_0} \right)^{b_2}\bigg|_{Upstream \ MTZ \ edge}+\frac{b_3}{1.01-\left( \frac{C}{C_0} \right)^{b_4}\bigg|_{Upstream \ MTZ \ edge}}`"
   "EBCT of the partial MTZ at the time of replacement", ":math:`EBCT_{MTZ} = \left( T_{MTZ}-T \right)EBCT_{min}`"
   "Length EBCT of the partial MTZ at the time of replacement", ":math:`EBCT_{MTZ}=\frac{L_{MTZ}}{v_s}`"
   "Saturation fraction of the bed at the time of replacement calculated by the trapezoid rule", ":math:`\bar{\left( \frac{q}{q_{e}} \right)}=\frac{\left( 1\left( L-L_{MTZ} \right) \right)+\frac{1}{2}\left( \frac{q}{q_{e}}\bigg|_{Upstream \ MTZ \ edge}+\frac{q}{q_{e}}\bigg|_{z=L}\right)\left( L_{MTZ} \right)}{L}`"
   "Length EBCT of the partial MTZ at the time of replacement", ":math:`EBCT_{MTZ}=\frac{L_{MTZ}}{v_s}`"

if ``film_transfer_coefficient_type`` or ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Hayduk-Laudie correlation", ":math:`D_l\left[ \frac{m^2}{s} \right] = \frac{13.26\times 10^{-9}}{\left( \mu_w\left[ cP \right] \right)^{1.14}\left( V_b \left[ \frac{cm^3}{mol} \right]\right)^{1.14}}`"

if ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Reynolds number for packed beds", ":math:`Re=\frac{\rho_w\phi d_pv_i}{\epsilon\mu_w}`"
   "Schmidt number for packed beds", ":math:`Sc=\frac{\mu_w}{\rho_wD_l}`"
   "Gnielinski correlation", ":math:`k_f=\frac{\left[ 1+1.5\left( 1-\epsilon \right) \right]D_l}{d_p}\left( 2+0.644Re^{\frac{1}{2}}Sc^{\frac{1}{3}} \right)`"

if ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Surface diffusion coefficient correlation", ":math:`D_s=\left( S\!P\!D\!F\!R \right)\left( \frac{\epsilon_pC_0D_l}{\rho_aq_e\tau_p} \right)`"

Class Documentation
-------------------

* :mod:`watertap.unit_models.gac`

References
-----------
Hand, D. W., Crittenden, J. C., & Thacker, W. E. (1984). Simplified models for design of fixed-bed adsorption systems. Journal
of Environmental Engineering, 110(2), 440-456.

Crittenden, J., Rhodes, R., Hand, D., Howe, K., & Tchobanoglous, G. (2012). MWHs Water Treatment. Principles and Design. Editorial
John Wiley & Sons.