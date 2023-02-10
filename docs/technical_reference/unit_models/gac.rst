Granular Activated Carbon (GAC)
===============================
This is an empirical, performance-based granular activated carbon (GAC) model that works under the following criteria and assumptions:
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

The batch operation results of the CPS are converted to approximate steady-state results for intuitive use of the model for
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
(i.e., temperature, pressure, component flowrates) that should be fixed for the model to be fully specified.
In association with using the Freundlich adsorption isotherm and empirical model, the following 9 variables are almost always fixed
and may be derived from experimental data:

   * Freundlich isotherm :math:`k` parameter
   * Freundlich isotherm :math:`\frac{1}{n}` parameter
   * Stanton number equation parameters :math:`a_0` and :math:`a_1`
   * throughput ratio equation parameters :math:`b_0`, :math:`b_1`, :math:`b_2`, :math:`b_3` and :math:`b_4`

Additionally, the following 9 variables are traditionally fixed:

   * the target dimensionless concentration *or* the dimensionless concentration at the time of replacement
   * empty bed contact time
   * bed voidage
   * superficial velocity *or* bed length
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
It also includes 1 StateBlock (adsorbed) for the solute that is adsorbed into the GAC particles.
The material removed in the adsorbed state block is simulated as liquid phase solute but should be interpreted as solute that has adsorbed
into the solid phase GAC particles. The steady state mass removal and replacement rate of the GAC itself is provided as a unit model
variable and excluded from flowsheet material balances. Therefore, GAC is never included as a component in property package specifications.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'Solute', 'Background solutes/ions]*"

\*Adsorbed "Solute" provided in the ``target_species`` argument of the unit model.
\*"Background solutes/ions" are provided in the imported property model.

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
   "Equilibrium concentration of adsorbed phase with liquid phase", ":math:`q_e`", "equil_conc", "None", ":math:`\left( \text{kg}_\text{adsorbate}\text{/kg}_\text{adsorbent} \right)`"
   "Mass of contaminant adsorbed at the time of replacement", ":math:`M_{solute}`", "mass_adsorbed", "None", ":math:`\text{kg}`"
   "Mass of contaminant adsorbed if fully saturated", ":math:`M_{solute\text{,}e}`", "mass_adsorbed_saturated", "None", ":math:`\text{kg}`"
   "Adsorber bed void fraction", ":math:`\epsilon`", "bed_voidage", "None", ":math:`\text{dimensionless}`"
   "Adsorber bed volume", ":math:`V`", "bed_volume", "None", ":math:`\text{m}^3`"
   "Adsorber bed area", ":math:`A`", "bed_area", "None", ":math:`\text{m}^2`"
   "Adsorber bed length", ":math:`L`", "bed_length", "None", ":math:`\text{m}`"
   "Mass of fresh GAC in the bed", ":math:`M_{GAC}`", "bed_mass_gac", "None", ":math:`\text{kg}`"
   "Superficial velocity", ":math:`v_s`", "velocity_sup", "None", ":math:`\text{m/s}`"
   "Interstitial velocity", ":math:`v_i`", "velocity_int", "None", ":math:`\text{m/s}`"
   "GAC particle porosity", ":math:`\epsilon_p`", "particle_porosity", "None", ":math:`\text{dimensionless}`"
   "GAC particle apparent density", ":math:`\rho_a`", "particle_dens_app", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle bulk density", ":math:`\rho_b`", "particle_dens_bulk", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle solid density", ":math:`\rho_s`", "particle_dens_sol", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle diameter", ":math:`d_p`", "particle_dia", "None", ":math:`\text{m}`"
   "Average dimensionless concentration of the effluent in the operating duration", ":math:`\frac{\bar{C}}{C_{0}}\bigg{|}_{z=L}`", "conc_ratio_avg", "None", ":math:`\text{dimensionless}`"
   "Dimensionless concentration of the effluent at the time of replacement", ":math:`\frac{C}{C_{0}}\bigg{|}_{z=L,\,t=t_{op}}`", "conc_ratio_replace", "None", ":math:`\text{dimensionless}`"
   "Approximate saturation of the GAC in the bed at the time of replacement", ":math:`\frac{\bar{q}}{q_{e}}\bigg{|}_{t=t_{op}}`", "gac_saturation_replace", "None", ":math:`\text{dimensionless}`"
   "Empty bed contact time", ":math:`EBCT`", "ebct", "None", ":math:`\text{s}`"
   "Mass throughput ratio", ":math:`T`", "mass_throughput", "None", ":math:`\text{dimensionless}`"
   "Residence time", ":math:`\tau`", "res_time", "None", ":math:`\text{s}`"
   "Elapsed operation time", ":math:`t_{op}`", "elap_time", "None", ":math:`\text{s}`"
   "Bed volumes treated", ":math:`BVT`", "bed_volumes_treated", "None", ":math:`\text{dimensionless}`"
   "Steady state GAC replacement rate", ":math:`\dot{m}_{GAC}`", "gac_mass_replace_rate", "None", ":math:`\text{m/s}`"
   "Liquid phase film transfer coefficient", ":math:`k_f`", "kf", "None", ":math:`\text{m/s}`"
   "Surface diffusion coefficient", ":math:`D_s`", "ds", "None", ":math:`\text{m}^2\text{/s}`"
   "Solute distribution parameter", ":math:`D_g`", "dg", "None", ":math:`\text{dimensionless}`"
   "Biot number", ":math:`Bi`", "N_Bi", "None", ":math:`\text{dimensionless}`"
   "Minimum Stanton number for CPS", ":math:`St_{min}`", "min_N_St", "None", ":math:`\text{dimensionless}`"
   "Minimum empty bed contact time for CPS", ":math:`EBCT_{min}`", "min_ebct", "None", ":math:`\text{s}`"
   "Minimum residence time for CPS", ":math:`\tau_{min}`", "min_res_time", "None", ":math:`\text{s}`"
   "Minimum elapsed operation time for CPS", ":math:`t_{min}`", "min_elap_time", "None", ":math:`\text{s}`"
   "Mass throughput ratio for the upstream edge of the MTZ", ":math:`T\left( \frac{C}{C_{0}}\bigg{|}_{Upstream\,MTZ\,edge} \right)`", "mass_throughput_mtz_upstream", "None", ":math:`\text{dimensionless}`"
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
   "Surface-to-pore diffusion flux ratio", ":math:`S\!P\!D\!F\!R`", "spdfr", "None", ":math:`\text{dimensionless}`"

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
   "Elapsed operation time", ":math:`t_{op}=t_{min}+\left( \tau-\tau_{min} \right)\left( D_g+1 \right)`"
   "Density relation to bed voidage", ":math:`\epsilon=1-\frac{\rho_b}{\rho_a}`"
   "Density relation to particle porosity", ":math:`\epsilon_p=1-\frac{\rho_a}{\rho_s}`"
   "Steady state GAC replacement rate", ":math:`\dot{m}_{GAC}=\frac{M_{GAC}}{t_{op}}`"
   "Adsorber bed volume", ":math:`EBCT=\frac{V}{Q}`"
   "Mass of GAC in a fresh adsorber bed", ":math:`M_{GAC}=V\rho_b`"
   "Velocity relationship", ":math:`v_i=\frac{v_s}{\epsilon}`"
   "Adsorbed bed area", ":math:`A=\frac{Q}{v_s}`"
   "Adsorbed bed length", ":math:`EBCT=\frac{L}{v_s}`"
   "Bed volumes treated", ":math:`BVT=\frac{t_{op}\epsilon}{\tau}`"
   "Mass of solute adsorbed if the bed was fully saturated", ":math:`M_{solute\text{,}e}=q_eM_{GAC}`"
   "Saturation fraction of the bed at the time of replacement", ":math:`\frac{\bar{q}}{q_{e}}\bigg{|}_{t=t_{op}}=\frac{M_{solute}}{M_{solute\text{,}e}}`"
   "Mass throughput ratio of the upstream edge of the MTZ", ":math:`T\big{|}_{Upstream\,MTZ\,edge}=b_0+b_1\left( \frac{C}{C_0}\bigg|_{Upstream\,MTZ\,edge} \right)^{b_2}+\frac{b_3}{1.01-\left( \frac{C}{C_0}\Big|_{Upstream\,MTZ\,edge} \right)^{b_4}}`"
   "EBCT of the partial MTZ at the time of replacement", ":math:`EBCT_{MTZ} = \left( T\big{|}_{Upstream\,MTZ\,edge}-T \right)EBCT_{min}`"
   "Length EBCT of the partial MTZ at the time of replacement", ":math:`EBCT_{MTZ}=\frac{L_{MTZ}}{v_s}`"
   "Saturation fraction of the bed at the time of replacement calculated by the trapezoid rule", ":math:`\frac{\bar{q}}{q_{e}}\bigg{|}_{t=t_{op}}=\frac{1\left( L-L_{MTZ} \right)+\frac{1}{2}\left( \frac{q}{q_{e}}\bigg|_{Upstream \ MTZ \ edge}+\frac{q}{q_{e}}\bigg|_{z=L}\right)\left( L_{MTZ} \right)}{L}`"

if ``film_transfer_coefficient_type`` or ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Hayduk-Laudie correlation*", ":math:`D_l\left[ \frac{m^2}{s} \right] = \frac{13.26\times 10^{-9}}{\left( \mu_w\left[ cP \right] \right)^{1.14}\left( V_b \left[ \frac{cm^3}{mol} \right]\right)^{1.14}}`"

\*Subscript :math:`w` denotes liquid phase properties, here those of pure water is used considering trace solute concentrations.

if ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Reynolds number for packed beds*", ":math:`Re=\frac{\rho_w\phi d_pv_i}{\epsilon\mu_w}`"
   "Schmidt number for packed beds*", ":math:`Sc=\frac{\mu_w}{\rho_wD_l}`"
   "Gnielinski correlation", ":math:`k_f=\frac{\left[ 1+1.5\left( 1-\epsilon \right) \right]D_l}{d_p}\left( 2+0.644Re^{\frac{1}{2}}Sc^{\frac{1}{3}} \right)`"

\*Subscript :math:`w` denotes liquid phase properties, here those of pure water is used considering trace solute concentrations.

if ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Surface diffusion coefficient correlation", ":math:`D_s=\left( S\!P\!D\!F\!R \right)\left( \frac{\epsilon_pC_0D_l}{\rho_aq_e\tau_p} \right)`"

Costing Method
---------------

Costing Method Variables
+++++++++++++++++++++++++

The following parameters are constructed when applying the GAC costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Default Value", "Units"

   "Number of GAC contactors in operation in parallel", ":math:`N_{op}`", "gac_num_contactors_op", "1", ":math:`\text{dimensionless}`"
   "Number of off-line redundant GAC contactors in parallel", ":math:`N_{red}`", "gac_num_contactors_redundant", "1", ":math:`\text{dimensionless}`"
   "GAC contactor polynomial cost coefficient 0", ":math:`x_0`", "gac_contactor_cost_coeff_0", "10010.9", ":math:`$`"
   "GAC contactor polynomial cost coefficient 1", ":math:`x_1`", "gac_contactor_cost_coeff_1", "2204.95", ":math:`$/m^3`"
   "GAC contactor polynomial cost coefficient 2", ":math:`x_2`", "gac_contactor_cost_coeff_2", "-15.9378", ":math:`$/\left( m^3 \right)^2`"
   "GAC contactor polynomial cost coefficient 3", ":math:`x_3`", "gac_contactor_cost_coeff_3", "0.110592", ":math:`$/\left( m^3 \right)^3`"
   "Reference maximum value of GAC mass needed for initial charge where economy of scale no longer discounts the unit price", ":math:`M_{GAC}^{ref}`", "bed_mass_gac_max_ref", "18143.7", ":math:`kg`"
   "GAC adsorbent exponential cost pre-exponential coefficient", ":math:`y_0`", "gac_adsorbent_unit_cost_coeff", "4.58342", ":math:`$/kg`"
   "GAC adsorbent exponential cost parameter coefficient", ":math:`y_1`", "gac_adsorbent_unit_cost_exp_coeff ", "-1.25311e-5", ":math:`kg^{-1}`"
   "GAC other cost power law coefficient", ":math:`z_0`", "gac_other_cost_coeff", "16660.7", ":math:`$/\left( m^3 \right)^{z_1}`"
   "GAC other cost power law exponent", ":math:`z_1`", "gac_other_cost_exp", "0.552207", ":math:`\text{dimensionless}`"
   "Fraction of spent GAC adsorbent that can be regenerated for reuse", ":math:`f_{regen}`", "gac_regen_frac", "0.70", ":math:`\text{dimensionless}`"
   "Unit cost to regenerate spent GAC adsorbent by an offsite regeneration facility", ":math:`C_{regen}`", "gac_regen_unit_cost", "4.28352", ":math:`$/kg`"
   "Unit cost to makeup spent GAC adsorbent with fresh adsorbent", ":math:`C_{makeup}`", "gac_makeup_unit_cost", "4.58223", ":math:`$/kg`"


Costing GAC contactors is defaulted to purchasing 1 operational and 1 redundant contactor for alternating operation. For large systems this may be a poor
assumption considering vessel sizing and achieving pseudo-steady state.  The number of contactors input by the user should justify reasonable
(commercially available) dimensions of identical modular contactors in parallel. When costing several operational vessels, the area reported
in the unit model should be interpreted as the sum of the areas across all operating GAC contactors.

The following variables are constructed when applying the GAC costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units"

   "Unit contactor(s) capital cost", ":math:`C_{cap,bed}`", "contactor_cost", ":math:`$`"
   "GAC adsorbent cost per unit mass", ":math:`C_{carbon}`", "adsorbent_unit_cost", ":math:`$/kg`"
   "Unit adsorbent capital cost", ":math:`C_{cap,carbon}`", "adsorbent_cost", ":math:`$`"
   "Unit other process capital cost", ":math:`C_{cap,other}`", "other_process_cost", ":math:`$`"
   "Cost to regenerate spent GAC adsorbent by an offsite regeneration facility", ":math:`C_{op,regen}`", "gac_regen_cost", ":math:`$/year`"
   "Cost to makeup spent GAC adsorbent with fresh adsorbent", ":math:`C_{op,makeup}`", "gac_makeup_cost", ":math:`$/year`"

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

Operating costs are calculated as the cost to replace spent GAC adsorbent in the contactor beds. Energy for backwash and booster
pumps are considered negligible compared to the regeneration costs. Since the replacement adsorbent purchases are expected to be
purchased in bulk at smaller quantities than the initial charge, the cost of fresh GAC adsorbent for makeup has an independent
cost per unit mass variable, expected to be higher than the initial charge unit cost.

    .. math::

        & C_{op,tot} = C_{op,regen}+C_{op,makeup} \\\\
        & C_{op,regen} = f_{regen}C_{regen}\dot{m}_{GAC} \\\\
        & C_{op,makeup} = \left( 1-f_{regen} \right)C_{makeup}\dot{m}_{GAC}

Code Documentation
-------------------

* :mod:`watertap.unit_models.gac`
* :meth:`watertap.costing.watertap_costing_package.WaterTAPCostingData.cost_gac`

References
-----------
Hand, D. W., Crittenden, J. C., & Thacker, W. E. (1984). Simplified models for design of fixed-bed adsorption systems.
Journal of Environmental Engineering, 110(2), 440-456.

Crittenden, J., Rhodes, R., Hand, D., Howe, K., & Tchobanoglous, G. (2012). MWHs Water Treatment. Principles and Design.
EditorialJohn Wiley & Sons.

United States Environmental Protection Agency. (2021). Work Breakdown Structure-Based Cost Model for Granular Activated
Carbon Drinking Water Treatment.