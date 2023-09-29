Granular Activated Carbon (GAC)
===============================
This is an empirical, performance-based granular activated carbon (GAC) model that works under the following criteria and assumptions:
   * simulation of this unit model is only supported with the Multi-Component Aqueous Solution (MCAS) property package
   * supports a single liquid phase only
   * supports adsorption of a single solute species only while other species are considered inert
   * supports steady-state only
   * assumes isothermal conditions
   * model performance is independent of a gravity-fed or pressurized GAC unit, therefore assumes isobaric conditions

.. index::
   pair: watertap.unit_models.gac;gac

.. currentmodule:: watertap.unit_models.gac

Introduction
------------
The implemented model for estimating GAC performance is adapted from a simplified model originally presented in Hand, 1984 and
further elaborated in Crittenden, 2012. This formulation is denoted as the constant-pattern homogeneous surface diffusion model (CPHSDM).
As a GAC system is operated as a batch process, a mass transfer zone (MTZ) is formed in the bed where a concentration profile, or
breakthrough curve, develops as a function of the adsorption properties. This MTZ is bounded by saturated GAC upstream and
fresh GAC downstream. The CPHSDM is valid under the assumption that the shape of the MTZ is constant as it travels through the
bed and a constant pattern solution (CPS) may be determined. The CPS is calculated through a multistep procedure utilizing common
dimensionless groups applied in polynomial fits to determine performance. Therefore, coefficients used in the polynomial must be
derived from experimental data of the intended system to produce valid results. Coefficients for common compounds treated by GAC may be
found in both Hand, 1984 and Crittenden, 2012. The model is estimated to have within 10% error and therefore may be applied to bed lengths
shorter than the minimum length determined by the CPHSDM within the error threshold (in addition to being applicable to bed lengths greater than the minimum length determined by the CPHSDM).

The batch operation results of the CPS are converted to approximate steady-state results for intuitive use of the model
for flowsheet purposes. A visualization of the transformation is provided in Figure 1. For a traditional breakthrough
curve the CPHSDM method calculates the single point, single ``conc_ratio_replace`` and ``operational_time``, highlighted on the
breakthrough curve. This operational time is the amount of elapsed time after startup that the bed is refreshed with new
GAC adsorbent. Steady state concentration can be analogous to all of the effluent in this operational time being stored
as holdup, therefore the average concentration ratio is significantly less than concentration ratio at the time of
bed replacement, as many days pass before the start of the breakthrough. To approximate the average effluent
concentration in this time frame, the breakthrough curve is numerically integrated with the trapezoid rule. The curve is
discretized with respect to the concentration ratio instead of the (traditionally done 'x' variable) operational time
due to simplicity of solving the model equations.

.. figure:: ../../_static/unit_models/gac.png
    :width: 800
    :align: center

    Figure 1. Discretization of the breakthrough curve for the steady state approximation. The ratio of the shaded area
    to the area left of the vertical dotted line corresponding to ``operational_time`` is the calculated
    ``conc_ratio_avg``. Expected values of ``conc_ratio_avg`` are often less than 0.25 depending on the
    ``conc_ratio_replace`` setpoint.

Degrees of Freedom
------------------
In the default configuration of the GAC unit model there are 17 degrees of freedom in addition to the inlet state variables
(i.e., temperature, pressure, component flowrates) that should be fixed for the model to be fully specified.
In association with using the Freundlich adsorption isotherm and empirical model, the following 9 variables are almost always fixed
and may be derived from experimental data:

   * Freundlich isotherm :math:`k` parameter
   * Freundlich isotherm :math:`\frac{1}{n}` parameter
   * Stanton number equation parameters :math:`a_0` and :math:`a_1` (Hand, 1984), (Crittenden, 1987)
   * throughput ratio equation parameters :math:`b_0`, :math:`b_1`, :math:`b_2`, :math:`b_3` and :math:`b_4` (Hand, 1984)

Additionally, the following 9 variables are traditionally fixed:

   * particle apparent density
   * particle diameter
   * empty bed contact time
   * bed voidage *or* particle bulk density
   * superficial velocity *or* bed length
   * effluent to inlet concentration ratio at operational time *or* steady state approximation of average effluent to inlet concentration ratio in operational time by trapezoid rule *or* bed volumes treated
   * surface diffusion coefficient
   * liquid phase film transfer coefficient

When setting the configuration options to calculate the surface diffusion coefficient and liquid phase film transfer coefficient,
these respective variables are no longer specified and 4 newly introduced variables must be fixed. This excludes new variables
or parameters that may be required to be specified within the property package when called by the GAC model. This is a net result of 19 degrees of freedom.
Newly utilized variables that must be fixed include:

   * shape correction factor
   * particle porosity
   * tortuosity of the path that the adsorbate must take as compared to the radius
   * surface-to-pore diffusion flux ratio

Model Structure
------------------
The GAC model consists of 1 ControlVolume0DBlock (process_flow) for the process flow of the water treatment train.
The process flow includes 2 StateBlocks (inlet and outlet) which are used for mass and momentum balances.
It also includes 1 StateBlock (adsorbed) for the solute that is adsorbed into the GAC particles.
The material removed in the adsorbed state block is simulated as liquid phase solute but should be interpreted as solute that has adsorbed
into the solid phase GAC particles. The steady state mass removal and replacement rate of the GAC itself is provided as a unit model
variable and excluded from flowsheet material balances. Therefore, GAC is never included as a component in property package specifications.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "time", ":math:`t`", "[0]"
   "phases", ":math:`p`", "['Liq']"
   "components", ":math:`j`", "['H2O', target_species, background solutes]*"
   "species adsorbed", ":math:`\text{target_species}`", "[target_species]"
   "inert species", ":math:`\text{inert_species}`", "['H2O', target_species, background solutes] - [target_species]"
   "number of discretized operational time elements used for steady state approximation", ":math:`\text{ele_disc}`", "[0:elements_ss_approx]"
   "number of discretized trapezoidal area terms for steady state approximation", ":math:`\text{ele_index}`", "[1:elements_ss_approx]"

| \* target_species is provided in the ``target_species`` argument of the unit model and corresponds to the single solute which is adsorbed.
| \* ``inert_species`` are the difference in ``component_list - target_species``.

.. _GAC_variables:

Variables
----------
Supporting only a single solute, variables concerning the adsorption of a species are respective to the target species
although not explicitly indexed. All other species are inert from a mass balance perspective, but effects of background
species to the adsorption of the target species may be modeled by adjusting the Freundlich isotherm parameters and other
variables in the model.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Freundlich isotherm k parameter", ":math:`k`", "freund_k", "None", ":math:`\left(\text{m}^3\text{/kg}\right)^\left( \frac{1}{n} \right)`"
   "Freundlich isotherm 1/n parameter", ":math:`\frac{1}{n}`", "freund_ninv", "None", ":math:`\text{dimensionless}`"
   "surface diffusion coefficient", ":math:`D_s`", "ds", "None", ":math:`\text{m}^2\text{/s}`"
   "liquid phase film transfer coefficient", ":math:`k_f`", "kf", "None", ":math:`\text{m/s}`"
   "equilibrium concentration of adsorbed phase with liquid phase", ":math:`q_e`", "equil_conc", "None", ":math:`\left( \text{kg}_\text{adsorbate}\text{/kg}_\text{adsorbent} \right)`"
   "solute distribution parameter", ":math:`D_g`", "dg", "None", ":math:`\text{dimensionless}`"
   "Biot number", ":math:`Bi`", "N_Bi", "None", ":math:`\text{dimensionless}`"
   "superficial velocity", ":math:`u_s`", "velocity_sup", "None", ":math:`\text{m/s}`"
   "interstitial velocity", ":math:`u_i`", "velocity_int", "None", ":math:`\text{m/s}`"
   "bed void fraction", ":math:`\epsilon`", "bed_voidage", "None", ":math:`\text{dimensionless}`"
   "bed length", ":math:`L`", "bed_length", "None", ":math:`\text{m}`"
   "bed diameter", ":math:`D`", "bed_diameter", "None", ":math:`\text{m}`"
   "bed area", ":math:`A`", "bed_area", "None", ":math:`\text{m}^2`"
   "bed volume", ":math:`V`", "bed_volume", "None", ":math:`\text{m}^3`"
   "empty bed contact time", ":math:`EBCT`", "ebct", "None", ":math:`\text{s}`"
   "fluid residence time in the bed", ":math:`\tau`", "residence_time", "None", ":math:`\text{s}`"
   "mass of fresh gac in the bed", ":math:`M_{GAC}`", "bed_mass_gac", "None", ":math:`\text{kg}`"
   "gac apparent density", ":math:`\rho_a`", "particle_dens_app", "None", ":math:`\text{kg/}\text{m}^3`"
   "gac bulk density", ":math:`\rho_b`", "particle_dens_bulk", "None", ":math:`\text{kg/}\text{m}^3`"
   "gac particle diameter", ":math:`d_p`", "particle_dia", "None", ":math:`\text{m}`"
   "Stanton equation parameter 0", ":math:`a_0`", "a0", "None", ":math:`\text{dimensionless}`"
   "Stanton equation parameter 1", ":math:`a_1`", "a1", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 0", ":math:`b_0`", "b0", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 1", ":math:`b_1`", "b1", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 2", ":math:`b_2`", "b2", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 3", ":math:`b_3`", "b3", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 4", ":math:`b_4`", "b4", "None", ":math:`\text{dimensionless}`"
   "minimum Stanton number to achieve a constant pattern solution", ":math:`St_{min}`", "min_N_St", "None", ":math:`\text{dimensionless}`"
   "minimum empty bed contact time to achieve a constant pattern solution", ":math:`EBCT_{min}`", "min_ebct", "None", ":math:`\text{s}`"
   "specific throughput from empirical equation", ":math:`T`", "throughput", "None", ":math:`\text{dimensionless}`"
   "minimum fluid residence time in the bed to achieve a constant pattern solution", ":math:`\tau_{min}`", "min_residence_time", "None", ":math:`\text{s}`"
   "minimum operational time of the bed from fresh to achieve a constant pattern solution", ":math:`t_{min}`", "min_operational_time", "None", ":math:`\text{s}`"
   "effluent to inlet concentration ratio at operational time", ":math:`\frac{C}{C_{0}}\bigg{|}_{z=L,/,t=t_{op}}`", "conc_ratio_replace", "None", ":math:`\text{dimensionless}`"
   "operational time of the bed from fresh", ":math:`t_{op}`", "operational_time", "None", ":math:`\text{s}`"
   "bed volumes treated at operational time", ":math:`BVT`", "bed_volumes_treated", "None", ":math:`\text{dimensionless}`"
   "specific throughput from empirical equation by discrete element", ":math:`T_{ele}`", "ele_throughput", "None", ":math:`x`"
   "minimum operational time of the bed from fresh to achieve a constant pattern solution by discrete element", ":math:`t_{min, ele}`", "ele_min_operational_time", "ele_index", ":math:`\text{s}`"
   "effluent to inlet concentration ratio at operational time by discrete element", ":math:`\left(\frac{C}{C_{0}}\right)_{ele}\bigg{|}_{z=L,/,t=t_{op_ e}}`", "ele_conc_ratio_replace", "ele_index", ":math:`\text{dimensionless}`"
   "operational time of the bed from fresh by discrete element", ":math:`t_{op, ele}`", "ele_operational_time", "ele_disc", ":math:`\text{s}`"
   "trapezoid rule of elements for numerical integration of average concentration ratio", ":math:`term_{ele}`", "ele_conc_ratio_avg", "ele_disc", ":math:`\text{dimensionless}`"
   "steady state approximation of average effluent to inlet concentration ratio in operational time by trapezoid rule", ":math:`\left(\frac{C}{C_{0}}\right)_{avg}`", "conc_ratio_avg", "None", ":math:`\text{dimensionless}`"
   "total mass of adsorbed species at operational time", ":math:`M`", "mass_adsorbed", "None", ":math:`\text{kg}`"
   "gac usage/replacement/regeneration rate", ":math:`\dot{m}_{GAC}`", "gac_usage_rate", "None", ":math:`\text{m/s}`"

The following variables are only built when specific configuration options are selected.

if ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Reynolds number", ":math:`Re`", "N_Re", "None", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc`", "N_Sc", "None", ":math:`\text{dimensionless}`"
   "shape correction factor", ":math:`SCF`", "shape_correction_factor", "None", ":math:`\text{dimensionless}`"

if ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "gac particle porosity", ":math:`\epsilon_p`", "particle_porosity", "None", ":math:`\text{dimensionless}`"
   "tortuosity of the path that the adsorbate must take as compared to the radius", ":math:`\tau_p`", "tort", "None", ":math:`\text{dimensionless}`"
   "surface-to-pore diffusion flux ratio", ":math:`S\!P\!D\!F\!R`", "spdfr", "None", ":math:`\text{dimensionless}`"

.. _GAC_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "equilibrium concentration", ":math:`q_e = kC_0^{1/n}`"
   "solute distribution parameter", ":math:`D_g=\frac{\rho_aq_e\left( 1-\epsilon \right)}{\epsilon C_0}`"
   "Biot number", ":math:`Bi=\frac{k_fd_p\left( 1-\epsilon \right)}{2D_sD_g\epsilon}`"
   "bed void fraction based on gac particle densities", ":math:`\epsilon=1-\frac{\rho_b}{\rho_a}`"
   "relating velocities based on bed voidage", ":math:`u_i=\frac{u_s}{\epsilon}`"
   "bed length based on velocity and ebct", ":math:`L=(EBCT)u_s`"
   "bed diameter and area relation", ":math:`A=\pi\left(\frac{D}{2}\right)^2`"
   "bed area based on velocity and volumetric flow", ":math:`A=\frac{Q}{u_s}`"
   "bed volume based on cylindrical dimensions", ":math:`V=AL`"
   "fluid residence time in the bed", ":math:`\tau=(EBCT)\epsilon`"
   "total mass of gac in the bed", ":math:`M_{GAC}=V\rho_b`"
   "minimum Stanton number to achieve constant pattern solution", ":math:`St_{min}=a_0Bi+a_1`"
   "minimum empty bed contact time to achieve constant pattern solution", ":math:`EBCT_{min}=\frac{St_{min}d_p}{2k_f\left( 1-\epsilon \right)}`"
   "throughput based on empirical 5-parameter regression", ":math:`T=b_0+b_1\left( \frac{C}{C_0} \right)^{b_2}+\frac{b_3}{1.01-\left( \frac{C}{C_0} \right)^{b_4}}`"
   "minimum fluid residence time in the bed to achieve a constant pattern solution", ":math:`\tau_{min}=EBCT_{min}\epsilon`"
   "minimum operational time of the bed from fresh to achieve a constant pattern solution", ":math:`t_{min}=\tau_{min}\left( D_g+1 \right)T`"
   "elapsed operational time between a fresh bed and the theoretical bed replacement", ":math:`t_{op}=t_{min}+\left( \tau-\tau_{min} \right)\left( D_g+1 \right)`"
   "bed volumes treated", ":math:`BVT=\frac{t_{op}\epsilon}{\tau}`"
   "throughput based on empirical 5-parameter regression by discretized element", ":math:`T_{ele}=b_0+b_1\left(\frac{C}{C_{0}}\right)_{ele}^{b_2}+\frac{b_3}{1.01-\left(\frac{C}{C_{0}}\right)_{ele}^{b_4}}`"
   "minimum operational time of the bed from fresh to achieve a constant pattern solution by discretized element", ":math:`t_{min, ele}=\tau_{min}\left( D_g+1 \right)T`"
   "creating evenly spaced discretized elements", ":math:`\frac{C}{C_{0}}\bigg{|}_{t=t_{op, ele}}=0.01+(ele-1)*\frac{\left(\frac{C}{C_{0}}\bigg{|}_{t=t_{op}}-0.01\right)}{num\text{_}ele}`"
   "finite element discretization of concentration ratios over time", ":math:`term_{ele}=\left(\frac{t_{op, ele}-t_{op, (ele-1)}}{t_{op}}\right)\frac{\left(\frac{C}{C_{0}}\bigg{|}_{t=t_{op, ele}}+{C_{0}}\bigg{|}_{t=t_{op, (ele-1)}}\right)}{2}`"
   "summation of finite elements for average concentration during operating time", ":math:`\left(\frac{C}{C_{0}}\right)_{avg}=\sum_{ele\text{_}index}term_{ele}`"
   "mass adsorbed in the operational time", ":math:`M=\frac{\dot{m}_{j}}{t_{op}}`"
   "steady state rate of new gac mass required", ":math:`\dot{m}_{GAC}=\frac{M_{GAC}}{t_{op}}`"

if ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Reynolds number calculation*", ":math:`Re=\frac{\rho_ld_pu_i}{\mu_l}`"
   "Schmidt number calculation*", ":math:`Sc=\frac{\mu_l}{\rho_sD_l}`"
   "liquid phase film transfer rate from the Gnielinshi correlation*", ":math:`k_f=(SCF)\frac{\left[ 1+1.5\left( 1-\epsilon \right) \right]D_l}{d_p}\left( 2+0.644Re^{\frac{1}{2}}Sc^{\frac{1}{3}} \right)`"

\*Subscript :math:`l` denotes bulk liquid phase properties, here those are supplied by the property package.

if ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "surface diffusion parameter (Crittenden, 1987)", ":math:`D_s=\left( S\!P\!D\!F\!R \right)\left( \frac{\epsilon_pC_0D_l}{\rho_aq_e\tau_p} \right)`"

Costing Method
---------------

Costing Method Variables
+++++++++++++++++++++++++

The following parameters are constructed when applying the GAC costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Default Value", "Units"

   "Number of GAC contactors in operation in parallel", ":math:`N_{op}`", "num_contactors_op", "1", ":math:`\text{dimensionless}`"
   "Number of off-line redundant GAC contactors in parallel", ":math:`N_{red}`", "num_contactors_redundant", "1", ":math:`\text{dimensionless}`"
   "Fraction of spent GAC adsorbent that can be regenerated for reuse", ":math:`f_{regen}`", "regen_frac", "0.70", ":math:`\text{dimensionless}`"
   "Reference maximum value of GAC mass needed for initial charge where economy of scale no longer discounts the unit price  (U.S. EPA, 2021)", ":math:`M_{GAC}^{ref}`", "bed_mass_gac_max_ref", "18143.7", ":math:`kg`"
   "Contactor polynomial cost coefficients", ":math:`x`", "contactor_cost_coeff", "tabulated", ":math:`\text{dimensionless}`"
   "Adsorbent exponential cost coefficients", ":math:`y`", "adsorbent_unit_cost_coeff", "tabulated", ":math:`\text{dimensionless}`"
   "Other process costs power law coefficients", ":math:`z`", "other_cost_param", "tabulated", ":math:`\text{dimensionless}`"
   "Unit cost to regenerate spent GAC adsorbent by an offsite regeneration facility", ":math:`C_{regen}`", "regen_unit_cost", "4.28352", ":math:`$/kg`"
   "Unit cost to makeup spent GAC adsorbent with fresh adsorbent", ":math:`C_{makeup}`", "makeup_unit_cost", "4.58223", ":math:`$/kg`"
   "Energy consumption polynomial coefficients", ":math:`alpha`", "energy_consumption_coeff", "tabulated", ":math:`\text{dimensionless}`"

Costing methods are available for steel pressure vessel contactors (default) and concrete gravity basin contactors. Given that the form of the costing component equations are different (polynomial, exponential, and power law), the units associated with the parameters are embedded in the constraints and not directly applied to the variable. Additionally, the index is generalized to its position ``([0:len(parameter_data)])`` in the list, although some parameters are coefficients while others are exponents (see equations below for details). Variables with the (U.S. EPA, 2021) citation are directly taken from previously determined expressions. Other variables are regressed from higher detailed costing methods in (U.S. EPA, 2021). The variations in costing parameters are tabulated below:

.. csv-table::
   :header: "Variable Name", "Contactor Type", "Index 0", "Index 1", "Index 2", "Index 3"

   "adsorbent_unit_cost_coeff (U.S. EPA, 2021)", "n/a", "4.58342", "-1.25311e-5", "", ""
   "contactor_cost_coeff (U.S. EPA, 2021)", "pressure", "10010.9", "2204.95", "-15.9378", "0.110592"
   "contactor_cost_coeff", "gravity", "75131.3", "735.550", "-1.01827", "0"
   "other_cost_param", "pressure", "16660.7", "0.552207", "", ""
   "other_cost_param", "gravity", "38846.9", "0.490571", "", ""
   "energy_consumption_coeff_data", "pressure", "8.09926e-4", "8.70577e-4", "0", ""
   "energy_consumption_coeff_data", "gravity", "0.123782", "0.132403", "-1.41512e-5", ""

Costing GAC contactors is defaulted to purchasing 1 operational and 1 redundant contactor for alternating operation. For large systems this may be a poor
assumption considering vessel sizing and achieving pseudo-steady state.  The number of contactors input by the user should justify reasonable
(commercially available) dimensions of identical modular contactors in parallel. When costing several operational vessels, the area reported
in the unit model should be interpreted as the sum of the areas across all operating GAC contactors. The costing
parameters may be selected from either steel pressure-fed vessels or concrete gravity-fed basins by the
``contactor_type`` argument. Note this only affects costing calculations. Volume dimensions calculations
within the model remain assuming a cylindrical bed.

The following variables are constructed when applying the GAC costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units"

   "Unit contactor(s) capital cost", ":math:`C_{cap,bed}`", "contactor_cost", ":math:`$`"
   "GAC adsorbent cost per unit mass", ":math:`C_{carbon}`", "adsorbent_unit_cost", ":math:`$/kg`"
   "Unit adsorbent capital cost", ":math:`C_{cap,carbon}`", "adsorbent_cost", ":math:`$`"
   "Unit other process capital cost", ":math:`C_{cap,other}`", "other_process_cost", ":math:`$`"
   "Cost to regenerate spent GAC adsorbent by an offsite regeneration facility", ":math:`C_{op,regen}`", "gac_regen_cost", ":math:`$/year`"
   "Cost to makeup spent GAC adsorbent with fresh adsorbent", ":math:`C_{op,makeup}`", "gac_makeup_cost", ":math:`$/year`"
   "Approximate GAC system energy consumption*", ":math:`P`", "energy_consumption", ":math:`kW`"

\*Energy consumption is the sum of energy required to operate booster, backwash, and residual pumps.

Capital Cost Calculations
+++++++++++++++++++++++++

Capital costs are determined by the summation of three costing terms. Each term is is calculated by a one parameter
(different for each term) function considering economy of scale.

    .. math::

        C_{cap,tot} = C_{cap,bed}+C_{cap,carbon}+C_{cap,other}

Contactor and GAC adsorbent capital costs are estimated using functions and parameters reported in US EPA, 2021. Contactors
are assumed to be carbon steel pressure vessels with plastic internals and are determined as a polynomial function of
individual contactor volume. The unit cost per kilogram of GAC adsorbent needed is calculated using an exponential
function. A maximum reference mass is imposed in the costing method to define a best available price where above
this required charge, the price would no longer be discounted. Other process costs (vessels, pipes, instrumentation,
and controls) included in the US EPA, 2021 model are aggregated into a separate term. The parameters for the power law
function with respect to the total system contactor volume were regressed using results from the US EPA, 2021 model.

    .. math::

        & C_{cap,bed} = \left( N_{op}+N_{red} \right)\left( x_0+x_1\left( \frac{V}{N_{op}} \right)+x_2\left( \frac{V}{N_{op}} \right)^2+x_3\left( \frac{V}{N_{op}} \right)^3 \right) \\\\
        & C_{carbon} = y_0e^{y_1M_{GAC}^{ref}} \\\\
        & C_{cap,carbon} = C_{carbon}M_{GAC} \\\\
        & C_{cap,other} = z_0\left( \left( N_{op}+N_{red} \right)\frac{V}{N_{op}} \right)^{z_1}

Note that given the the ability to alter the parameters in these correlations, GAC adsorbent unit costs (:math:`C_{carbon}`)
may be fixed to a value (:math:`y_0`) by setting :math:`y_1=0`.

Operating Cost Calculations
+++++++++++++++++++++++++++

Operating costs are calculated as the cost to replace spent GAC adsorbent in the contactor beds. Energy is costed as a
flow term by the WaterTAP costing method. Since the replacement adsorbent purchases are expected to be purchased in bulk
at smaller quantities than the initial charge, the cost of fresh GAC adsorbent for makeup has an independent cost per
unit mass variable, expected to be higher than the initial charge unit cost.

    .. math::

        & C_{op,tot} = C_{op,regen}+C_{op,makeup} \\\\
        & C_{op,regen} = f_{regen}C_{regen}\dot{m}_{GAC} \\\\
        & C_{op,makeup} = \left( 1-f_{regen} \right)C_{makeup}\dot{m}_{GAC} \\\\
        & P = \alpha_0+\alpha_1V+\alpha_2V^2

Code Documentation
-------------------

* :mod:`watertap.unit_models.gac`
* :mod:`watertap.costing.unit_models.gac`

References
-----------
Hand, D. W., Crittenden, J. C., & Thacker, W. E. (1984). Simplified models for design of fixed-bed adsorption systems.
Journal of Environmental Engineering, 110(2), 440-456.

Crittenden, J., Rhodes, R., Hand, D., Howe, K., & Tchobanoglous, G. (2012). MWHs Water Treatment. Principles and Design.
John Wiley & Sons.

Crittenden, J. C., Berrigan, J. K., Hand, D. W., & Lykins, B. (1987). Design of Rapid Fixed‐Bed Adsorption Tests for
Nonconstant Diffusivities. Journal of Environmental Engineering, 113(2), 243–259.

United States Environmental Protection Agency. (2021). Work Breakdown Structure-Based Cost Model for Granular Activated
Carbon Drinking Water Treatment. https://www.epa.gov/system/files/documents/2022-03/gac-documentation-.pdf_0.pdf